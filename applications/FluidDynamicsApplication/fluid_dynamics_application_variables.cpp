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

#include "fluid_dynamics_application_variables.h"
#include "custom_utilities/statistics_data.h"
#include "custom_utilities/statistics_record.h"

namespace Kratos
{
KRATOS_CREATE_VARIABLE(int,PATCH_INDEX)
KRATOS_CREATE_VARIABLE(double,TAUONE)
KRATOS_CREATE_VARIABLE(double,TAUTWO)
KRATOS_CREATE_VARIABLE(double,PRESSURE_MASSMATRIX_COEFFICIENT)
KRATOS_CREATE_VARIABLE(Vector, FLUID_STRESS)
KRATOS_CREATE_VARIABLE(Vector, GAPS)
KRATOS_CREATE_VARIABLE(double, DIVERGENCE)
KRATOS_CREATE_VARIABLE(double, FS_PRESSURE_GRADIENT_RELAXATION_FACTOR)

KRATOS_CREATE_VARIABLE(double,SUBSCALE_PRESSURE)
KRATOS_CREATE_VARIABLE(double, C_DES)

KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(SUBSCALE_VELOCITY)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(COARSE_VELOCITY)

KRATOS_CREATE_VARIABLE(double,FIC_BETA)
KRATOS_CREATE_VARIABLE(double,VOLUME_ERROR)
KRATOS_CREATE_VARIABLE(double, CONVECTION_SCALAR)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(CONVECTION_VELOCITY)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(CONVECTION_SCALAR_GRADIENT)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(AUXILIAR_VECTOR_VELOCITY)

// Darcy's flow variables
KRATOS_CREATE_VARIABLE(double, RESISTANCE)

// Wall modelling
KRATOS_CREATE_VARIABLE(bool, SLIP_TANGENTIAL_CORRECTION_SWITCH)

// Outlet inflow contribution
KRATOS_CREATE_VARIABLE(double, CHARACTERISTIC_VELOCITY)
KRATOS_CREATE_VARIABLE(bool, OUTLET_INFLOW_CONTRIBUTION_SWITCH)

// Adjoint variables
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(ADJOINT_FLUID_VECTOR_1)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(ADJOINT_FLUID_VECTOR_2)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(ADJOINT_FLUID_VECTOR_3)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(AUX_ADJOINT_FLUID_VECTOR_1)
KRATOS_CREATE_VARIABLE(double, ADJOINT_FLUID_SCALAR_1)
KRATOS_CREATE_VARIABLE(Vector, PRIMAL_RELAXED_SECOND_DERIVATIVE_VALUES)

// Non-Newtonian constitutive relations
KRATOS_CREATE_VARIABLE(double, REGULARIZATION_COEFFICIENT)

KRATOS_CREATE_VARIABLE(double, BINGHAM_SMOOTHER)
KRATOS_CREATE_VARIABLE(double, GEL_STRENGTH )

// Q-Criterion (for vortex visualization)
KRATOS_CREATE_VARIABLE(double, Q_VALUE)

// Vorticity
KRATOS_CREATE_VARIABLE(double, VORTICITY_MAGNITUDE)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(RECOVERED_PRESSURE_GRADIENT)

// For swimming DEM
KRATOS_CREATE_VARIABLE(Vector, NODAL_WEIGHTS)

// Embedded fluid variables
KRATOS_CREATE_VARIABLE(int, EMBEDDED_IS_ACTIVE)
KRATOS_CREATE_VARIABLE(double, SLIP_LENGTH)
KRATOS_CREATE_VARIABLE(double, EMBEDDED_WET_PRESSURE)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(EMBEDDED_WET_VELOCITY)

// Compressible fluid variables
KRATOS_CREATE_VARIABLE(bool, SHOCK_CAPTURING_SWITCH)
KRATOS_CREATE_VARIABLE(double, MASS_SOURCE)
KRATOS_CREATE_VARIABLE(double, HEAT_SOURCE)
KRATOS_CREATE_VARIABLE(double, HEAT_CAPACITY_RATIO)
KRATOS_CREATE_VARIABLE(double, REACTION_DENSITY)
KRATOS_CREATE_VARIABLE(double, REACTION_ENERGY)
KRATOS_CREATE_VARIABLE(double, MACH)
KRATOS_CREATE_VARIABLE(double, SHOCK_SENSOR)
KRATOS_CREATE_VARIABLE(double, SHEAR_SENSOR)
KRATOS_CREATE_VARIABLE(double, THERMAL_SENSOR)
KRATOS_CREATE_VARIABLE(double, ARTIFICIAL_MASS_DIFFUSIVITY)
KRATOS_CREATE_VARIABLE(double, ARTIFICIAL_CONDUCTIVITY)
KRATOS_CREATE_VARIABLE(double, ARTIFICIAL_BULK_VISCOSITY)
KRATOS_CREATE_VARIABLE(double, ARTIFICIAL_DYNAMIC_VISCOSITY)
KRATOS_CREATE_VARIABLE(double, DENSITY_PROJECTION)
KRATOS_CREATE_VARIABLE(double, TOTAL_ENERGY_PROJECTION)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DENSITY_GRADIENT)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(MOMENTUM_PROJECTION)
KRATOS_CREATE_VARIABLE(double, NUMERICAL_ENTROPY)

// Turbulence statistics
KRATOS_CREATE_VARIABLE( StatisticsRecord::Pointer, STATISTICS_CONTAINER)
KRATOS_CREATE_VARIABLE( StatisticsData, TURBULENCE_STATISTICS_DATA)
KRATOS_CREATE_VARIABLE( double, UPDATE_STATISTICS )

// Auxiliary variables
KRATOS_CREATE_VARIABLE(Matrix, VELOCITY_GRADIENT)
KRATOS_CREATE_VARIABLE(double, VELOCITY_DIVERGENCE)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(VELOCITY_ROTATIONAL)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DRAG_FORCE_CENTER)
KRATOS_CREATE_VARIABLE( double, SMOOTHING_COEFFICIENT )
KRATOS_CREATE_VARIABLE(double,SPECTRAL_RADIUS_LIMIT)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( INLET_NORMAL )

// Two-phase flow with surface tension
KRATOS_CREATE_VARIABLE( double, SURFACE_TENSION_COEFFICIENT )
KRATOS_CREATE_VARIABLE( bool, SURFACE_TENSION )
KRATOS_CREATE_VARIABLE( double, CURVATURE )

// Derivative variables
KRATOS_CREATE_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_COMPONENTS( STRAIN_RATE_2D )
KRATOS_CREATE_SYMMETRIC_3D_TENSOR_VARIABLE_WITH_COMPONENTS( STRAIN_RATE_3D )
KRATOS_CREATE_3D_TENSOR_VARIABLE_WITH_COMPONENTS( VELOCITY_GRADIENT_TENSOR )
KRATOS_CREATE_3D_TENSOR_VARIABLE_WITH_COMPONENTS( SHAPE_SENSITIVITY_GRADIENT_TENSOR )

// Two-phase flow momentum correction
KRATOS_CREATE_VARIABLE( bool, MOMENTUM_CORRECTION )
KRATOS_CREATE_VARIABLE( double, DISTANCE_CORRECTION )

KRATOS_CREATE_VARIABLE( double, Y_PLUS )

}
