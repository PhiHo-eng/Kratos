## Dictionary of the constant parameters used in the Variational Formulation
import sympy
from dataclasses import dataclass
from KratosMultiphysics.FluidDynamicsApplication.symbolic_generation.compressible_navier_stokes.src.defines \
    import CompressibleNavierStokesDefines as defs

@dataclass
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

@dataclass
class ShockCapturingParameters:
    def __init__(self):
        self.alpha = sympy.Symbol('data.alpha_sc', positive = True) # Artificial density diffusivity for shock capturing
        self.mu    = sympy.Symbol('data.mu_sc', positive = True) # Artificial dynamic viscosity for shock capturing
        self.beta  = sympy.Symbol('data.beta_sc', positive = True) # Artificial bulk viscosity for shock capturing
        self.lamb  = sympy.Symbol('data.lamb_sc', positive = True) # Artificial thermal conductivity for shock capturing

@dataclass
class ShockCapturingNodalParameters:
    def __init__(self, geometry):
        self.alpha = defs.Vector('data.alpha_sc_nodes', geometry.nnodes, positive=True) # Nodal artificial mass diffusivity
        self.mu    = defs.Vector('data.mu_sc_nodes', geometry.nnodes, positive=True)    # Nodal artificial dynamic viscosity
        self.beta  = defs.Vector('data.beta_sc_nodes', geometry.nnodes, positive=True)  # Nodal artificial bulk viscosity
        self.lamb  = defs.Vector('data.lamb_sc_nodes', geometry.nnodes, positive=True)  # Nodal artificial bulk viscosity

class PrimitiveMagnitudes:
    def __init__(self, geometry):
        self.P = sympy.Symbol('pressure', real=True)
        self.T = sympy.Symbol('temperature', real=True)
        self.V = defs.Vector('velocity', geometry.ndims, real=True)
        self.nnodes = geometry.nnodes
        self.ndims = geometry.ndims

    def Interpolate(self, mode, U_nodes, U_gauss, N_gauss, params):
        if mode == "gaussian":
            return self._GaussianInterpolation(U_gauss, params)
        elif mode == "nodal":
            return self._NodalInterpolation(U_nodes, N_gauss, params)

    def _GaussianInterpolation(self, U_gauss, params):
        mom_g = sympy.Matrix(U_gauss[1:-1])
        rho_g = U_gauss[0]
        e_tot_g = U_gauss[-1]

        V_g = self._velocity(rho_g, mom_g)
        T_g = self._temperature(rho_g, V_g, e_tot_g, params)
        P_g = self._pressure(rho_g, T_g, params)

        return (P_g, V_g, T_g)

    def _NodalInterpolation(self, U_nodes, N_gauss, params):
        V_nodes = defs.ZeroMatrix(self.nnodes, self.ndims)
        T_nodes = defs.ZeroVector(self.nnodes)
        P_nodes = defs.ZeroVector(self.nnodes)

        for n in range(self.nnodes):
            mom_n = U_nodes[n, 1:-1]
            rho_n = U_nodes[n, 0]
            e_tot_n = U_nodes[n, -1]

            V_nodes[n, :] = self._velocity(rho_n, mom_n)
            T_nodes[n] = self._temperature(rho_n, V_nodes[n, :], e_tot_n, params)
            P_nodes[n] = self._pressure(rho_n, T_nodes[n], params)

        V_g = V_nodes.transpose() * N_gauss
        T_g = T_nodes.transpose() * N_gauss
        P_g = P_nodes.transpose() * N_gauss

        V_g.simplify()
        T_g.simplify()
        P_g.simplify()

        # Asserts to avoid unintelligible errors from the depths of sympy
        assert(sympy.shape(V_g) == (self.ndims,1))
        assert(sympy.shape(T_g) == (1,1))
        assert(sympy.shape(P_g) == (1,1))

        return (P_g[0, 0], V_g, T_g[0,0])

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
