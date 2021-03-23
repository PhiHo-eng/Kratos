//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//
//


#include "particle_mechanics_application_variables.h"

namespace Kratos
{
    // Unused variables
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( NODAL_INTERNAL_FORCE )

    // Element
    KRATOS_CREATE_VARIABLE( int, MP_MATERIAL_ID )
    KRATOS_CREATE_VARIABLE( int, PARTICLES_PER_ELEMENT )
    KRATOS_CREATE_VARIABLE( int, MP_SUB_POINTS)
    KRATOS_CREATE_VARIABLE( double, MP_MASS )
    KRATOS_CREATE_VARIABLE( double, AUX_MASS)
    KRATOS_CREATE_VARIABLE( double, MP_DENSITY )
    KRATOS_CREATE_VARIABLE( double, MP_VOLUME )
    KRATOS_CREATE_VARIABLE( double, MP_RADIUS)
    KRATOS_CREATE_VARIABLE( double, MP_POTENTIAL_ENERGY )
    KRATOS_CREATE_VARIABLE( double, MP_KINETIC_ENERGY )
    KRATOS_CREATE_VARIABLE( double, MP_STRAIN_ENERGY )
    KRATOS_CREATE_VARIABLE( double, MP_TOTAL_ENERGY )
    KRATOS_CREATE_VARIABLE( double, MP_PRESSURE )
    KRATOS_CREATE_VARIABLE( double, PRESSURE_REACTION )
    KRATOS_CREATE_VARIABLE( double, MP_EQUIVALENT_STRESS)
    KRATOS_CREATE_VARIABLE( double, MP_DELTA_PLASTIC_STRAIN )
    KRATOS_CREATE_VARIABLE( double, MP_DELTA_PLASTIC_VOLUMETRIC_STRAIN )
    KRATOS_CREATE_VARIABLE( double, MP_DELTA_PLASTIC_DEVIATORIC_STRAIN )
    KRATOS_CREATE_VARIABLE( double, MP_EQUIVALENT_PLASTIC_STRAIN )
    KRATOS_CREATE_VARIABLE( double, MP_EQUIVALENT_PLASTIC_STRAIN_RATE)
    KRATOS_CREATE_VARIABLE( double, MP_ACCUMULATED_PLASTIC_VOLUMETRIC_STRAIN )
    KRATOS_CREATE_VARIABLE( double, MP_ACCUMULATED_PLASTIC_DEVIATORIC_STRAIN )
    KRATOS_CREATE_VARIABLE( double, MP_DAMAGE )
    KRATOS_CREATE_VARIABLE( double, MP_TEMPERATURE)
    KRATOS_CREATE_VARIABLE( double, NODAL_MPRESSURE )
    KRATOS_CREATE_VARIABLE(bool, IS_COMPRESSIBLE)

    // Constitutive Law
    KRATOS_CREATE_VARIABLE( ConstitutiveLaw::Pointer, CONSTITUTIVE_LAW_POINTER )
    // CL: Solid
    KRATOS_CREATE_VARIABLE( double, RAYLEIGH_ALPHA )
    KRATOS_CREATE_VARIABLE( double, RAYLEIGH_BETA )
    // CL: Mohr Coulomb
    KRATOS_CREATE_VARIABLE( double, COHESION )
    KRATOS_CREATE_VARIABLE( double, INTERNAL_DILATANCY_ANGLE )
    // CL: Mohr Coulomb Strain Softening
    KRATOS_CREATE_VARIABLE( double, INTERNAL_FRICTION_ANGLE_RESIDUAL )
    KRATOS_CREATE_VARIABLE( double, COHESION_RESIDUAL )
    KRATOS_CREATE_VARIABLE( double, INTERNAL_DILATANCY_ANGLE_RESIDUAL )
    KRATOS_CREATE_VARIABLE( double, SHAPE_FUNCTION_BETA )
    // CL: Johnson Cook
    KRATOS_CREATE_VARIABLE( double, REFERENCE_STRAIN_RATE)
    KRATOS_CREATE_VARIABLE( double, TAYLOR_QUINNEY_COEFFICIENT)
    KRATOS_CREATE_VARIABLE( double, MP_HARDENING_RATIO)
    KRATOS_CREATE_VARIABLE( double, FRACTURE_TOUGHNESS)
    KRATOS_CREATE_VARIABLE( double, JC_D1)
    KRATOS_CREATE_VARIABLE(double, JC_D2)
    KRATOS_CREATE_VARIABLE(double, JC_D3)
    KRATOS_CREATE_VARIABLE(double, JC_D4)
    KRATOS_CREATE_VARIABLE(double, JC_D5)
    // CL: d+/d- Damage
    KRATOS_CREATE_VARIABLE(double, DAMAGE_TENSION)
    KRATOS_CREATE_VARIABLE(double, UNIAXIAL_STRESS_TENSION)
    KRATOS_CREATE_VARIABLE(double, UNIAXIAL_STRAIN_TENSION)
    KRATOS_CREATE_VARIABLE(double, UNIAXIAL_DAMAGED_STRESS_TENSION)
    KRATOS_CREATE_VARIABLE(double, THRESHOLD_TENSION)
    KRATOS_CREATE_VARIABLE(double, DAMAGE_COMPRESSION)
    KRATOS_CREATE_VARIABLE(double, UNIAXIAL_STRESS_COMPRESSION)
    KRATOS_CREATE_VARIABLE(double, UNIAXIAL_STRAIN_COMPRESSION)
    KRATOS_CREATE_VARIABLE(double, UNIAXIAL_DAMAGED_STRESS_COMPRESSION)
    KRATOS_CREATE_VARIABLE(double, THRESHOLD_COMPRESSION)
    KRATOS_CREATE_VARIABLE(double, YIELD_STRESS_TENSION)
    KRATOS_CREATE_VARIABLE(double, FRACTURE_ENERGY_TENSION)
    KRATOS_CREATE_VARIABLE(double, DAMAGE_ONSET_STRESS_COMPRESSION)
    KRATOS_CREATE_VARIABLE(double, YIELD_STRESS_COMPRESSION)
    KRATOS_CREATE_VARIABLE(double, RESIDUAL_STRESS_COMPRESSION)
    KRATOS_CREATE_VARIABLE(double, YIELD_STRAIN_COMPRESSION)
    KRATOS_CREATE_VARIABLE(double, FRACTURE_ENERGY_COMPRESSION)
    KRATOS_CREATE_VARIABLE(double, BIAXIAL_COMPRESSION_MULTIPLIER)
    KRATOS_CREATE_VARIABLE(double, INTEGRATION_IMPLEX)
    KRATOS_CREATE_VARIABLE(double, BEZIER_CONTROLLER_S1)
    KRATOS_CREATE_VARIABLE(double, BEZIER_CONTROLLER_EP1)
    KRATOS_CREATE_VARIABLE(double, BEZIER_CONTROLLER_EP2)
    KRATOS_CREATE_VARIABLE(double, BEZIER_CONTROLLER_EP3)
    KRATOS_CREATE_VARIABLE(double, BEZIER_CONTROLLER_EP4)
    KRATOS_CREATE_VARIABLE(double, SHEAR_COMPRESSION_REDUCTOR)
    KRATOS_CREATE_VARIABLE(double, STRAIN_RATE_FACTOR_C1_TENSION)
    KRATOS_CREATE_VARIABLE(double, STRAIN_RATE_FACTOR_C2_TENSION)
    KRATOS_CREATE_VARIABLE(double, STRAIN_RATE_FACTOR_C1_COMPRESSION)
    KRATOS_CREATE_VARIABLE(double, STRAIN_RATE_FACTOR_C2_COMPRESSION)
    KRATOS_CREATE_VARIABLE(double, STRAIN_RATE_FACTOR_C1_YOUNGS_MOD)
    KRATOS_CREATE_VARIABLE(double, STRAIN_RATE_FACTOR_C2_YOUNGS_MOD)
    KRATOS_CREATE_VARIABLE(double, PLASTICITY_FACTOR_B_MINUS)
    KRATOS_CREATE_VARIABLE(int, TENSION_YIELD_MODEL)

    // Mesh variables
    KRATOS_CREATE_VARIABLE(std::vector<typename Geometry<Node<3>>::Pointer>, GEOMETRY_NEIGHBOURS)

    // Conditions
    // Essential Boundary Conditions
    KRATOS_CREATE_VARIABLE( double, PENALTY_FACTOR )
    KRATOS_CREATE_VARIABLE( int, MPC_BOUNDARY_CONDITION_TYPE )

    // Nodal load variables
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(POINT_LOAD)
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(LINE_LOAD)
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(SURFACE_LOAD)

    // MP element variable
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( MP_COORD )
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( MP_DISPLACEMENT )
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( MP_VELOCITY )
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( MP_ACCELERATION )
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( MP_VOLUME_ACCELERATION )
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(AUX_MOMENTA)
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(AUX_INERTIA)
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(AUX_RESIDUAL)
    KRATOS_CREATE_VARIABLE( Vector, MP_CAUCHY_STRESS_VECTOR )
    KRATOS_CREATE_VARIABLE( Vector, MP_ALMANSI_STRAIN_VECTOR )

    // MP condition variable
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( MPC_COORD )
    KRATOS_CREATE_VARIABLE( int, MPC_CONDITION_ID )
    KRATOS_CREATE_VARIABLE( bool, MPC_IS_NEUMANN )
    KRATOS_CREATE_VARIABLE( double, MPC_AREA )
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( MPC_NORMAL )
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( MPC_DISPLACEMENT )
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( MPC_IMPOSED_DISPLACEMENT )
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( MPC_VELOCITY )
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( MPC_IMPOSED_VELOCITY )
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( MPC_ACCELERATION )
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( MPC_IMPOSED_ACCELERATION )
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( MPC_CONTACT_FORCE )
    KRATOS_CREATE_VARIABLE( int, PARTICLES_PER_CONDITION )

    // Grid node variable
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( NODAL_MOMENTUM)
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( NODAL_INERTIA)

    // Solver related variables
    KRATOS_CREATE_VARIABLE(bool, IGNORE_GEOMETRIC_STIFFNESS)
    KRATOS_CREATE_VARIABLE(bool, IS_AXISYMMETRIC)
    KRATOS_CREATE_VARIABLE(bool, IS_PQMPM)
    KRATOS_CREATE_VARIABLE(bool, IS_MAKE_NORMAL_MP_IF_PQMPM_FAILS)
    KRATOS_CREATE_VARIABLE(double, PQMPM_SUBPOINT_MIN_VOLUME_FRACTION)

    // Explicit time integration variables
    KRATOS_CREATE_VARIABLE(bool, CALCULATE_MUSL_VELOCITY_FIELD)
    KRATOS_CREATE_VARIABLE(bool, IS_EXPLICIT)
    KRATOS_CREATE_VARIABLE(bool, IS_EXPLICIT_CENTRAL_DIFFERENCE)
    KRATOS_CREATE_VARIABLE(int, EXPLICIT_STRESS_UPDATE_OPTION)
    KRATOS_CREATE_VARIABLE(bool, CALCULATE_EXPLICIT_MP_STRESS)
    KRATOS_CREATE_VARIABLE(bool, EXPLICIT_MAP_GRID_TO_MP)
    KRATOS_CREATE_VARIABLE(bool, IS_FIX_EXPLICIT_MP_ON_GRID_EDGE)
    KRATOS_CREATE_VARIABLE(bool, EXPLICIT_CONTACT_RELEASE)
    KRATOS_CREATE_VARIABLE(std::string, EXPLICIT_CONTACT_RELEASE_MODEL_PART)

    // Cosim coupling variables
    KRATOS_CREATE_VARIABLE(int, INTERFACE_EQUATION_ID)
    KRATOS_CREATE_VARIABLE(bool, IS_COSIM_COUPLED)
}