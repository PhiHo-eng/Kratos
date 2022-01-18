## Dictionary of the constant parameters used in the Variational Formulation
import sympy
import KratosMultiphysics.sympy_fe_utilities as KratosSympy
from KratosMultiphysics.FluidDynamicsApplication.symbolic_generation.compressible_navier_stokes.src.defines \
    import CompressibleNavierStokesDefines as defs

class FormulationParameters:
    def __init__(self, geometry):
        self.mu      = sympy.Symbol('data.mu', positive = True)       # Dynamic viscosity
        self.h       = sympy.Symbol('data.h', positive = True)        # Element size
        self.lamb    = sympy.Symbol('data.lambda', positive = True)   # Thermal Conductivity of the fluid
        self.c_v     = sympy.Symbol('data.c_v', positive = True)      # Specific Heat at Constant volume
        self.gamma   = sympy.Symbol('data.gamma',positive = True)     # Gamma (Cp/Cv)
        self.stab_c1 = sympy.Symbol('stab_c1', positive = True)  # Algorithm constant
        self.stab_c2 = sympy.Symbol('stab_c2', positive = True)  # Algorithm constant
        self.stab_c3 = sympy.Symbol('stab_c3', positive = True)  # Algorithm constant
        self.dim = geometry.ndims

class ShockCapturingParameters:
    def __init__(self):
        self.alpha = sympy.Symbol('data.alpha_sc', positive = True) # Artificial density diffusivity for shock capturing
        self.mu    = sympy.Symbol('data.mu_sc', positive = True) # Artificial dynamic viscosity for shock capturing
        self.beta  = sympy.Symbol('data.beta_sc', positive = True) # Artificial bulk viscosity for shock capturing
        self.lamb  = sympy.Symbol('data.lamb_sc', positive = True) # Artificial thermal conductivity for shock capturing

class ShockCapturingNodalParameters:
    def __init__(self, geometry):
        self.alpha = defs.Vector('data.alpha_sc_nodes', geometry.nnodes, positive=True) # Nodal artificial mass diffusivity
        self.mu    = defs.Vector('data.mu_sc_nodes', geometry.nnodes, positive=True)    # Nodal artificial dynamic viscosity
        self.beta  = defs.Vector('data.beta_sc_nodes', geometry.nnodes, positive=True)  # Nodal artificial bulk viscosity
        self.lamb  = defs.Vector('data.lamb_sc_nodes', geometry.nnodes, positive=True)  # Nodal artificial bulk viscosity

class PrimitiveMagnitudes:
    """
    Primitive variables and their gradients.
    Gradients defined as:
    ```
    grad_f[i, j] := df_j/dx_i.
    ```

    There are two interpolation methods.

    For a given `f:U->V`, where `U` is the
    set of conservative variables and `V` is the set of primitive ones,

    1. Nodal interpolation computes V at the nodes and then interpolates it:
    ```
    V(x) = Σ_i N_i(x)·f(U(x_i))     for i=1,2...nnodes
    ```

    2. Gaussian interpolation interpolates U and then computes V at the gauss points:
    ```
    V(x) = f(Σ_i N_i(x)·U(x_i))     for i=1,2...nnodes
    ```
    """

    GAUSSIAN = 1
    NODAL = 2

    def __init__(self, geometry, params, Ug, grad_Ug, mode):
        if mode == "nodal":
            self.mode = self.NODAL
            self.P = sympy.Symbol('pressure', real=True)
            self.T = sympy.Symbol('temperature', real=True)
            self.V = defs.Vector('velocity', geometry.ndims, real=True)
            self.grad_P = defs.Matrix('grad_pressure', 1, geometry.ndims, real=True)
            self.grad_V = defs.Matrix('grad_velocity', geometry.ndims, geometry.ndims, real=True)
            self.grad_T = defs.Matrix('grad_temperature', 1, geometry.ndims, real=True)
        elif mode == "gaussian":
            self.mode = self.GAUSSIAN
            (self.P, self.V, self.T) = self._PrimitivesFromConservatives(Ug, params)
            (self.grad_P, self.grad_V, self.grad_T) = self._PrimitiveGradientsFromConservatives(Ug, grad_Ug, params)
            self.nnodes = geometry.nnodes
            self.ndims = geometry.ndims
        else:
            raise ValueError("Unrecognized magnitude interpolation mode.")

    def InterpolateAndSubstitute(self, expr, U_nodes, N_gauss, DN_DX_gauss, params):
        """
        If necessary, the primitive variable symbols are substituted by their
        expressions at the gauss points in terms of the conservative variables.

        Arguments
        ---------
        expr:        The expression to perform the substitution in.
        U_nodes:     The values of the conservative variables at the nodes.
        N_gauss:     The values of the shape functions at this gauss point.
        DN_DX_gauss: The values of the shape functions' gradients at this gauss point.
        params:      The physical parameters of the formulation.
        """
        if self.mode == self.GAUSSIAN:
            # Interpolation done during __init__
            return

        V_nodes = defs.ZeroMatrix(self.nnodes, self.ndims)
        T_nodes = defs.ZeroVector(self.nnodes)
        P_nodes = defs.ZeroVector(self.nnodes)

        for n in range(self.nnodes):
            (P, V, T) = self._PrimitivesFromConservatives(U_nodes[n, :], params)
            P_nodes[n] = P
            V_nodes[n, :] = V
            T_nodes[n] = T

        P_g = P_nodes.transpose() * N_gauss
        V_g = V_nodes.transpose() * N_gauss
        T_g = T_nodes.transpose() * N_gauss

        # Asserts to avoid unintelligible errors from the depths of sympy
        assert(sympy.shape(V_g) == (self.ndims,1))
        assert(sympy.shape(T_g) == (1,1))
        assert(sympy.shape(P_g) == (1,1))

        P_g = P_g[0, 0]
        T_g = T_g[0, 0]

        # Gradients
        grad_P_g = KratosSympy.DfiDxj(DN_DX_gauss, P_nodes)
        grad_V_g = KratosSympy.DfiDxj(DN_DX_gauss, V_nodes)
        grad_T_g = KratosSympy.DfiDxj(DN_DX_gauss, T_nodes)

        # Substitution
        KratosSympy.SubstituteScalarValue(expr, self.P, P_g)
        KratosSympy.SubstituteMatrixValue(expr, self.V, V_g)
        KratosSympy.SubstituteScalarValue(expr, self.T, T_g)
        KratosSympy.SubstituteMatrixValue(expr, self.grad_P, grad_P_g)
        KratosSympy.SubstituteMatrixValue(expr, self.grad_V, grad_V_g)
        KratosSympy.SubstituteMatrixValue(expr, self.grad_T, grad_T_g)

    @classmethod
    def _PrimitivesFromConservatives(cls, U, params):
        rho = U[0]
        mom = sympy.Matrix(U[1:-1])
        e_tot = U[-1]

        V = cls._velocity(rho, mom)
        T = cls._temperature(rho, V, e_tot, params)
        P = cls._pressure(rho, T, params)

        return (P, V, T)

    @classmethod
    def _PrimitiveGradientsFromConservatives(cls, U, grad_U, params):
        rho = U[0]
        mom = sympy.Matrix(U[1:-1])
        e_tot = U[-1]

        V = cls._velocity(rho, mom)
        T = cls._temperature(rho, V, e_tot, params)

        grad_rho = sympy.Matrix(grad_U[0, :])                   # d rho / dx_i
        grad_mom = sympy.Matrix(grad_U[1:-1, :])                # d V_j / dx_i
        grad_e_tot = sympy.Matrix(grad_U[-1, :])                # d e_tot / dx_i

        grad_V = cls._velocity_gradient(rho, mom, grad_rho, grad_mom)
        grad_T = cls._temperature_gradient(rho, e_tot, grad_rho, grad_e_tot, V, grad_V, params)
        grad_P = cls._pressure_gradient(rho, grad_rho, T, grad_T, params)

        return (grad_P.T, grad_V.T, grad_T.T)

    @classmethod
    def _velocity(cls, rho, mom):
        return mom/rho

    @classmethod
    def _temperature(cls, rho, vel, e_tot, params):
        e_kinetic = 0.5 * rho * sum([v**2 for v in vel])
        return (e_tot - e_kinetic) / (rho * params.c_v)

    @classmethod
    def _pressure(cls, rho, T, params):
        R = (params.gamma - 1) * params.c_v
        return rho * R * T

    @classmethod
    def _velocity_gradient(cls, rho, mom, grad_rho, grad_mom):
        """
        Velocity gradient. Gradients defined as:

        grad_f := df_j/dx_i

        """
        return (grad_mom*rho - mom*grad_rho) / rho**2


    @classmethod
    def _temperature_gradient(cls, rho, e_tot, grad_rho, grad_e_tot, vel, grad_vel, params):
        """
        Temperature gradient.  Gradients defined as:

        grad_f := df_j/dx_i

        """
        gradient_kinetic_energy = vel.transpose() * grad_vel
        gradient_total_energy   = (rho*grad_e_tot - e_tot*grad_rho) / rho**2
        return (gradient_total_energy - gradient_kinetic_energy) / params.c_v

    @classmethod
    def _pressure_gradient(cls, rho, grad_rho, T, grad_T, params):
        """
        Pressure gradient.  Gradients defined as:

        grad_f := df_j/dx_i

        """
        R = (params.gamma - 1) * params.c_v
        return R * (grad_rho*T + rho*grad_T)
