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


#if !defined(KRATOS_FLUID_DYNAMICS_APPLICATION_VARIABLES_H_INCLUDED )
#define  KRATOS_FLUID_DYNAMICS_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"
#include "includes/dem_variables.h"  //TODO: must be removed eventually

#include "custom_utilities/statistics_record.h"
#include "custom_utilities/statistics_data.h"

namespace Kratos
{
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, int,PATCH_INDEX)
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double,TAUONE)
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double,TAUTWO)
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double,PRESSURE_MASSMATRIX_COEFFICIENT)
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double,SUBSCALE_PRESSURE)
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double,C_DES)
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double,CHARACTERISTIC_VELOCITY)
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( FLUID_DYNAMICS_APPLICATION, SUBSCALE_VELOCITY)
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( FLUID_DYNAMICS_APPLICATION, COARSE_VELOCITY)
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double,FIC_BETA)
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, Vector, FLUID_STRESS )
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, Vector, GAPS )
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double, DIVERGENCE)
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double, AUX_DISTANCE)
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double, FS_PRESSURE_GRADIENT_RELAXATION_FACTOR)

// Adjoint variables
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( FLUID_DYNAMICS_APPLICATION, ADJOINT_FLUID_VECTOR_1 )
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( FLUID_DYNAMICS_APPLICATION, ADJOINT_FLUID_VECTOR_2 )
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( FLUID_DYNAMICS_APPLICATION, ADJOINT_FLUID_VECTOR_3 )
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( FLUID_DYNAMICS_APPLICATION, AUX_ADJOINT_FLUID_VECTOR_1 )
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( FLUID_DYNAMICS_APPLICATION, RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUX_ADJOINT_FLUID_VECTOR_1 )
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double, ADJOINT_FLUID_SCALAR_1 )
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, Vector, PRIMAL_RELAXED_SECOND_DERIVATIVE_VALUES )

// Non-Newtonian constitutive relations
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double, REGULARIZATION_COEFFICIENT)

// To be removed
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double, BINGHAM_SMOOTHER)
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double, GEL_STRENGTH)

// Q-Criterion (for vortex visualization)
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double,Q_VALUE)

// Vorticity
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double,VORTICITY_MAGNITUDE)
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(FLUID_DYNAMICS_APPLICATION, RECOVERED_PRESSURE_GRADIENT)
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, Vector,NODAL_WEIGHTS)

// Embedded fluid variables
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, int, EMBEDDED_IS_ACTIVE)
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double, SLIP_LENGTH)
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double, PENALTY_COEFFICIENT)
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double, EMBEDDED_WET_PRESSURE)
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( FLUID_DYNAMICS_APPLICATION, EMBEDDED_WET_VELOCITY)

// Compressible fluid variable
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, bool, SHOCK_CAPTURING_SWITCH)
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double, MASS_SOURCE)
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double, HEAT_SOURCE)
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double, HEAT_CAPACITY_RATIO)
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double, REACTION_DENSITY)
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double, REACTION_ENERGY)
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double, MACH)
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double, SHOCK_SENSOR)
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double, SHEAR_SENSOR)
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double, THERMAL_SENSOR)
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double, ARTIFICIAL_CONDUCTIVITY)
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double, ARTIFICIAL_BULK_VISCOSITY)
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double, ARTIFICIAL_DYNAMIC_VISCOSITY)
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double, DENSITY_PROJECTION)
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double, TOTAL_ENERGY_PROJECTION)
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( FLUID_DYNAMICS_APPLICATION, DENSITY_GRADIENT)
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( FLUID_DYNAMICS_APPLICATION, MOMENTUM_PROJECTION)

// Turbulence statistics
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, StatisticsRecord::Pointer, STATISTICS_CONTAINER)
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, StatisticsData, TURBULENCE_STATISTICS_DATA)
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double, UPDATE_STATISTICS )

// Auxiliary variables
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, Matrix, VELOCITY_GRADIENT)
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double, VELOCITY_DIVERGENCE)
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( FLUID_DYNAMICS_APPLICATION, VELOCITY_ROTATIONAL)
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( FLUID_DYNAMICS_APPLICATION, DRAG_FORCE_CENTER )
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double, SMOOTHING_COEFFICIENT )

// Two-phase flow with surface tension
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double, SURFACE_TENSION_COEFFICIENT )
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, bool, SURFACE_TENSION )
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double, CURVATURE )

// Derivative variables
KRATOS_DEFINE_SYMMETRIC_2D_TENSOR_APPLICATION_VARIABLE_WITH_COMPONENTS( FLUID_DYNAMICS_APPLICATION, STRAIN_RATE_2D )
KRATOS_DEFINE_SYMMETRIC_3D_TENSOR_APPLICATION_VARIABLE_WITH_COMPONENTS( FLUID_DYNAMICS_APPLICATION, STRAIN_RATE_3D )
KRATOS_DEFINE_3D_TENSOR_APPLICATION_VARIABLE_WITH_COMPONENTS( FLUID_DYNAMICS_APPLICATION, VELOCITY_GRADIENT_TENSOR )
KRATOS_DEFINE_3D_TENSOR_APPLICATION_VARIABLE_WITH_COMPONENTS( FLUID_DYNAMICS_APPLICATION, SHAPE_SENSITIVITY_GRADIENT_TENSOR )

KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, Vector, RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUXILIARY_1 )
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, Vector, RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUXILIARY_2 )
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double, RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUXILIARY_1_EXTREMUM )
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double, RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUXILIARY_2_EXTREMUM )
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, bool, COMPUTE_RESPONSE_FUNCTION_INTERPOLATION_ERROR )

}

#endif	/* KRATOS_FLUID_DYNAMICS_APPLICATION_VARIABLES_H_INCLUDED */
