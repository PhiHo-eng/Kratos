//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define_python.h"
#include "fluid_dynamics_application.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_python/add_custom_constitutive_laws_to_python.h"
#include "custom_python/add_custom_response_functions_to_python.h"


namespace Kratos
{

namespace Python
{



PYBIND11_MODULE(KratosFluidDynamicsApplication,m)
{
    namespace py = pybind11;

    py::class_<KratosFluidDynamicsApplication,
           KratosFluidDynamicsApplication::Pointer,
           KratosApplication >(m,"KratosFluidDynamicsApplication")
           .def(py::init<>())
           ;

    AddCustomConstitutiveLawsToPython(m);
    AddCustomStrategiesToPython(m);
    AddCustomUtilitiesToPython(m);
    AddCustomProcessesToPython(m);
    AddCustomResponseFunctionsToPython(m);

    //registering variables in python
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,PATCH_INDEX);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,TAUONE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,TAUTWO);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,PRESSURE_MASSMATRIX_COEFFICIENT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,Y_WALL);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,SUBSCALE_PRESSURE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,C_DES);
//        KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,C_SMAGORINSKY);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,CHARACTERISTIC_VELOCITY);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,GAPS);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,DIVERGENCE);

    // For Non-Newtonian constitutive relations
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,REGULARIZATION_COEFFICIENT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,BINGHAM_SMOOTHER);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,GEL_STRENGTH);

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,SUBSCALE_VELOCITY);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,COARSE_VELOCITY);

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,FIC_BETA);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,VOLUME_ERROR);
    // Adjoint variables
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, ADJOINT_FLUID_VECTOR_1 )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, ADJOINT_FLUID_VECTOR_2 )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, ADJOINT_FLUID_VECTOR_3 )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, AUX_ADJOINT_FLUID_VECTOR_1 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ADJOINT_FLUID_SCALAR_1 )

    // Embedded fluid variables
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,EMBEDDED_IS_ACTIVE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,SLIP_LENGTH);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,PENALTY_COEFFICIENT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,EMBEDDED_WET_PRESSURE);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,EMBEDDED_WET_VELOCITY);

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,Q_VALUE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,VORTICITY_MAGNITUDE);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,RECOVERED_PRESSURE_GRADIENT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,NODAL_WEIGHTS);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,AUX_DISTANCE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,FS_PRESSURE_GRADIENT_RELAXATION_FACTOR)

    // Compressible fluid variables
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,SHOCK_CAPTURING_SWITCH);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,MASS_SOURCE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,HEAT_SOURCE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,HEAT_CAPACITY_RATIO);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,REACTION_DENSITY);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,REACTION_ENERGY);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,MACH);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,SHOCK_SENSOR);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,SHEAR_SENSOR);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,THERMAL_SENSOR);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,ARTIFICIAL_MASS_DIFFUSIVITY);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,ARTIFICIAL_CONDUCTIVITY);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,ARTIFICIAL_BULK_VISCOSITY);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,ARTIFICIAL_DYNAMIC_VISCOSITY);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,DENSITY_GRADIENT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,DENSITY_PROJECTION);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,TOTAL_ENERGY_PROJECTION);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,MOMENTUM_PROJECTION);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NUMERICAL_ENTROPY);

    // Two-phase flow with surface tension
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,SURFACE_TENSION_COEFFICIENT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,SURFACE_TENSION);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,CURVATURE);

    // Two-phase flow momentum correction
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,MOMENTUM_CORRECTION )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,MOMENTUM_MASS_CORRECTION )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,DISTANCE_CORRECTION);

    // Auxiliary variables
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,VELOCITY_GRADIENT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, VELOCITY_DIVERGENCE);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, VELOCITY_ROTATIONAL);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,SMOOTHING_COEFFICIENT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,SPECTRAL_RADIUS_LIMIT)
}


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
