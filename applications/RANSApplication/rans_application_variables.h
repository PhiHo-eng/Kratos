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

#if !defined(KRATOS_RANS_APPLICATION_VARIABLES_H_INCLUDED)
#define KRATOS_RANS_APPLICATION_VARIABLES_H_INCLUDED

// System includes
#include <vector>
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "containers/variable.h"

namespace Kratos
{
    // incompressible potential flow specific variables
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, double, VELOCITY_POTENTIAL )
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, double, PRESSURE_POTENTIAL )
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, int, RANS_IS_INLET )
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, int, RANS_IS_OUTLET )
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, int, RANS_IS_STRUCTURE )
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, int, RANS_IS_STEADY )

    // residual based flux corrected stabilization variables
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, double, RANS_STABILIZATION_DISCRETE_UPWIND_OPERATOR_COEFFICIENT )
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, double, RANS_STABILIZATION_DIAGONAL_POSITIVITY_PRESERVING_COEFFICIENT )

    // algebraic flux corrected stabilization variables
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, double, AFC_POSITIVE_ANTI_DIFFUSIVE_FLUX )
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, double, AFC_NEGATIVE_ANTI_DIFFUSIVE_FLUX )
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, double, AFC_POSITIVE_ANTI_DIFFUSIVE_FLUX_LIMIT )
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, double, AFC_NEGATIVE_ANTI_DIFFUSIVE_FLUX_LIMIT )

    // k-epsilon-high-re turbulence modelling variables
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, double, TURBULENT_KINETIC_ENERGY )
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, double, TURBULENT_ENERGY_DISSIPATION_RATE )
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, double, TURBULENT_KINETIC_ENERGY_RATE )
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, double, TURBULENT_ENERGY_DISSIPATION_RATE_2 )
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, double, RANS_AUXILIARY_VARIABLE_1 )
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, double, RANS_AUXILIARY_VARIABLE_2 )
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, double, TURBULENCE_RANS_C_MU )
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, double, TURBULENT_KINETIC_ENERGY_SIGMA )
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, double, TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA )
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, double, TURBULENCE_RANS_C1 )
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, double, TURBULENCE_RANS_C2 )

    // k-omega turbulence modelling specific additional variables
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, double, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE )
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, double, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_2 )
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, double, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA )
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, double, TURBULENCE_RANS_BETA )
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, double, TURBULENCE_RANS_GAMMA )

    // k-omega-sst turbulence modelling specific additional variables
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, double, TURBULENT_KINETIC_ENERGY_SIGMA_1 )
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, double, TURBULENT_KINETIC_ENERGY_SIGMA_2 )
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, double, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA_1 )
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, double, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA_2 )
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, double, TURBULENCE_RANS_A1 )
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, double, TURBULENCE_RANS_BETA_1 )
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, double, TURBULENCE_RANS_BETA_2 )
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, double, VON_KARMAN )

    // wall function condition specific additional variables
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, double, RANS_Y_PLUS )
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, Vector, GAUSS_RANS_Y_PLUS )
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, double, RANS_LINEAR_LOG_LAW_Y_PLUS_LIMIT )
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, double, WALL_SMOOTHNESS_BETA )
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, double, WALL_CORRECTION_FACTOR )
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, int, RANS_IS_WALL_FUNCTION_ACTIVE )
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( RANS_APPLICATION, FRICTION_VELOCITY )

    // formulation specific variables
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, std::vector<std::string>, ANALYSIS_STEPS )
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, std::string, WALL_MODEL_PART_NAME )
    KRATOS_DEFINE_APPLICATION_VARIABLE( RANS_APPLICATION, double, NUMBER_OF_NEIGHBOUR_CONDITIONS )

    // adjoint variables
    KRATOS_DEFINE_APPLICATION_VARIABLE(RANS_APPLICATION, double, RANS_SCALAR_1_ADJOINT_1 )
    KRATOS_DEFINE_APPLICATION_VARIABLE(RANS_APPLICATION, double, RANS_SCALAR_1_ADJOINT_2 )
    KRATOS_DEFINE_APPLICATION_VARIABLE(RANS_APPLICATION, double, RANS_SCALAR_1_ADJOINT_3 )
    KRATOS_DEFINE_APPLICATION_VARIABLE(RANS_APPLICATION, double, RANS_AUX_ADJOINT_SCALAR_1 )

    KRATOS_DEFINE_APPLICATION_VARIABLE(RANS_APPLICATION, double, RANS_SCALAR_2_ADJOINT_1 )
    KRATOS_DEFINE_APPLICATION_VARIABLE(RANS_APPLICATION, double, RANS_SCALAR_2_ADJOINT_2 )
    KRATOS_DEFINE_APPLICATION_VARIABLE(RANS_APPLICATION, double, RANS_SCALAR_2_ADJOINT_3 )
    KRATOS_DEFINE_APPLICATION_VARIABLE(RANS_APPLICATION, double, RANS_AUX_ADJOINT_SCALAR_2 )

    // primal solution location storage variables
    KRATOS_DEFINE_APPLICATION_VARIABLE(RANS_APPLICATION, std::string, RANS_PRIMAL_SOLUTION_LOCATION_1 )

    // gauss point quantity variables
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(RANS_APPLICATION, RANS_GAUSS_EFFECTIVE_VELOCITY )
    KRATOS_DEFINE_APPLICATION_VARIABLE(RANS_APPLICATION, double, RANS_GAUSS_EFFECTIVE_KINEMATIC_VISCOSITY )
    KRATOS_DEFINE_APPLICATION_VARIABLE(RANS_APPLICATION, double, RANS_GAUSS_REACTION_TERM )
    KRATOS_DEFINE_APPLICATION_VARIABLE(RANS_APPLICATION, double, RANS_GAUSS_SOURCE_TERM )

} // namespace Kratos

#endif /* KRATOS_RANS_APPLICATION_VARIABLES_H_INCLUDED */
