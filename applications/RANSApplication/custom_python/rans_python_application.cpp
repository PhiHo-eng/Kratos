//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//


// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define.h"
#include "rans_application.h"
#include "rans_application_variables.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"


namespace Kratos {
namespace Python {

PYBIND11_MODULE(KratosRANSApplication,m)
{
    namespace py = pybind11;

    py::class_<KratosRANSApplication,
        KratosRANSApplication::Pointer,
        KratosApplication>(m, "KratosRANSApplication")
        .def(py::init<>())
        ;

    AddCustomStrategiesToPython(m);
    AddCustomUtilitiesToPython(m);
    AddCustomProcessesToPython(m);

    // incompressible potential flow specific variables
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, VELOCITY_POTENTIAL )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PRESSURE_POTENTIAL )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RANS_IS_INLET )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RANS_IS_OUTLET )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RANS_IS_STRUCTURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RANS_IS_STEADY )

    // residual based flux corrected stabilization variables
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RANS_STABILIZATION_DISCRETE_UPWIND_OPERATOR_COEFFICIENT )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RANS_STABILIZATION_DIAGONAL_POSITIVITY_PRESERVING_COEFFICIENT )

    // algebraic flux corrected stabilization variables
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, AFC_POSITIVE_ANTI_DIFFUSIVE_FLUX )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, AFC_NEGATIVE_ANTI_DIFFUSIVE_FLUX )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, AFC_POSITIVE_ANTI_DIFFUSIVE_FLUX_LIMIT )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, AFC_NEGATIVE_ANTI_DIFFUSIVE_FLUX_LIMIT )

    // k-epsilon-high-re turbulence modelling variables
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TURBULENT_KINETIC_ENERGY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TURBULENT_ENERGY_DISSIPATION_RATE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TURBULENT_KINETIC_ENERGY_RATE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TURBULENT_ENERGY_DISSIPATION_RATE_2 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RANS_AUXILIARY_VARIABLE_1 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RANS_AUXILIARY_VARIABLE_2 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TURBULENT_KINETIC_ENERGY_SIGMA )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TURBULENCE_RANS_C_MU )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TURBULENCE_RANS_C1 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TURBULENCE_RANS_C2 )

    // k-omega turbulence modelling specific additional variables
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_2 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TURBULENCE_RANS_BETA )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TURBULENCE_RANS_GAMMA )

    // k-omega-sst turbulence modelling specific additional variables
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TURBULENT_KINETIC_ENERGY_SIGMA_1 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TURBULENT_KINETIC_ENERGY_SIGMA_2 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA_1 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA_2 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TURBULENCE_RANS_A1 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TURBULENCE_RANS_BETA_1 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TURBULENCE_RANS_BETA_2 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, VON_KARMAN )

    // wall function condition specific additional variables
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RANS_Y_PLUS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RANS_LINEAR_LOG_LAW_Y_PLUS_LIMIT )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, WALL_SMOOTHNESS_BETA )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, WALL_CORRECTION_FACTOR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RANS_IS_WALL_FUNCTION_ACTIVE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, FRICTION_VELOCITY )

    // formulation specific variables
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, WALL_MODEL_PART_NAME )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NUMBER_OF_NEIGHBOUR_CONDITIONS )

    // adjoint variables
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RANS_SCALAR_1_ADJOINT_1 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RANS_SCALAR_1_ADJOINT_2 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RANS_SCALAR_1_ADJOINT_3 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RANS_AUX_ADJOINT_SCALAR_1 )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RANS_SCALAR_2_ADJOINT_1 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RANS_SCALAR_2_ADJOINT_2 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RANS_SCALAR_2_ADJOINT_3 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RANS_AUX_ADJOINT_SCALAR_2 )

    // primal solution location storage variables
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RANS_PRIMAL_SOLUTION_LOCATION_1)

    // gauss point quantity variables
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RANS_GAUSS_EFFECTIVE_VELOCITY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RANS_GAUSS_EFFECTIVE_KINEMATIC_VISCOSITY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RANS_GAUSS_REACTION_TERM )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RANS_GAUSS_SOURCE_TERM )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RANS_GAUSS_STABILIZATION_TAU )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RANS_GAUSS_SCALAR_CONSISTENCY_MULTIPLIER )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RANS_GAUSS_K_OMEGA_SST_F1 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RANS_GAUSS_K_OMEGA_SST_BLENDED_GAMMA )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RANS_GAUSS_K_OMEGA_SST_F1_ARG_1 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RANS_GAUSS_K_OMEGA_SST_F1_ARG_2 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RANS_GAUSS_K_OMEGA_SST_F1_ARG_3 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RANS_GAUSS_K_OMEGA_SST_F1_ARG_4 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RANS_GAUSS_K_OMEGA_SST_F1_ARG_5 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RANS_GAUSS_K_OMEGA_SST_F1_ARG_6 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RANS_GAUSS_PRODUCTION_TERM )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RANS_GAUSS_VELOCITY_GRADIENT )
}

} // namespace Python.
} // namespace Kratos.

#endif // KRATOS_PYTHON defined
