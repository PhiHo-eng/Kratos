// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//


// System includes


// External includes

// Project includes
#include "includes/define.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/quadrilateral_2d_8.h"
#include "geometries/quadrilateral_2d_9.h"
#include "geometries/quadrilateral_interface_2d_4.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/tetrahedra_3d_10.h"
#include "geometries/prism_interface_3d_6.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_20.h"
#include "geometries/hexahedra_3d_27.h"
#include "geometries/hexahedra_interface_3d_8.h"
#include "geometries/point_2d.h"
#include "geometries/point_3d.h"
#include "geometries/line_2d_2.h"
#include "geometries/line_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/triangle_3d_6.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/quadrilateral_3d_8.h"
#include "geometries/quadrilateral_3d_9.h"
#include "geometries/quadrilateral_interface_3d_4.h"
#include "includes/variables.h"


// Application includes
#include "geo_mechanics_application.h"
#include "geo_mechanics_application_variables.h"


namespace Kratos {
// We define the node type
typedef Node<3> NodeType;

KratosGeoMechanicsApplication::KratosGeoMechanicsApplication():
    KratosApplication("GeoMechanicsApplication")

    // small strain elements

    // drained small strain elements

    // undrained small strain elements

    // small strain FIC elements

    // small strain different order elements

    // small strain interface elements

    // Updated-Lagrangian elements

    // Updated-Lagrangian different order elements

    // geo structural elements

    // conditions

    {}

void KratosGeoMechanicsApplication::Register() {
    // calling base class register to register Kratos components
    KratosApplication::Register();
    KRATOS_INFO("") << " KRATOS___                             \n"
                    << "     //   ) )                          \n"
                    << "    //         ___      ___            \n"
                    << "   //  ____  //___) ) //   ) )         \n"
                    << "  //    / / //       //   / /          \n"
                    << " ((____/ / ((____   ((___/ /  MECHANICS\n"
                    << " Initializing KratosGeoMechanicsApplication..." << std::endl;



    //Register Elements
    // Small strain elements

    // Drained small strain elements

    // Undrained small strain elements

    // Small strain FIC elements

    // Small strain different order elements

    // Small strain interface elements

    // Updated-Lagranian elements

    // Updated-Lagrangian different order elements

    // Register geo structural elements

    //Register Conditions

    //Register Constitutive Laws

    //Register Variables
    KRATOS_REGISTER_VARIABLE( VELOCITY_COEFFICIENT )
    KRATOS_REGISTER_VARIABLE( DT_PRESSURE_COEFFICIENT )

    KRATOS_REGISTER_VARIABLE( DT_WATER_PRESSURE )
    KRATOS_REGISTER_VARIABLE( NORMAL_FLUID_FLUX )

    KRATOS_REGISTER_VARIABLE( DENSITY_SOLID )
    KRATOS_REGISTER_VARIABLE( BULK_MODULUS_SOLID )
    KRATOS_REGISTER_VARIABLE( BULK_MODULUS_FLUID )

    KRATOS_REGISTER_VARIABLE( K0_MAIN_DIRECTION )
    KRATOS_REGISTER_VARIABLE( K0_VALUE_XX )
    KRATOS_REGISTER_VARIABLE( K0_VALUE_YY )
    KRATOS_REGISTER_VARIABLE( K0_VALUE_ZZ )

    KRATOS_REGISTER_VARIABLE( PERMEABILITY_XX )
    KRATOS_REGISTER_VARIABLE( PERMEABILITY_YY )
    KRATOS_REGISTER_VARIABLE( PERMEABILITY_ZZ )
    KRATOS_REGISTER_VARIABLE( PERMEABILITY_XY )
    KRATOS_REGISTER_VARIABLE( PERMEABILITY_YZ )
    KRATOS_REGISTER_VARIABLE( PERMEABILITY_ZX )

    KRATOS_REGISTER_VARIABLE( MINIMUM_JOINT_WIDTH )
    KRATOS_REGISTER_VARIABLE( TRANSVERSAL_PERMEABILITY )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( FLUID_FLUX_VECTOR )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( LOCAL_FLUID_FLUX_VECTOR )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( LOCAL_STRESS_VECTOR )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( LOCAL_RELATIVE_DISPLACEMENT_VECTOR )
    KRATOS_REGISTER_VARIABLE( PERMEABILITY_MATRIX )
    KRATOS_REGISTER_VARIABLE( LOCAL_PERMEABILITY_MATRIX )

    KRATOS_REGISTER_VARIABLE( CRITICAL_DISPLACEMENT )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( TOTAL_DISPLACEMENT )

    KRATOS_REGISTER_VARIABLE( IS_CONVERGED )

    KRATOS_REGISTER_VARIABLE( TOTAL_STRESS_TENSOR )
    KRATOS_REGISTER_VARIABLE( TOTAL_STRESS_VECTOR )

    KRATOS_REGISTER_VARIABLE( CAUCHY_STRAIN_TENSOR )
    KRATOS_REGISTER_VARIABLE( CAUCHY_STRAIN_VECTOR )

    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE )
    KRATOS_REGISTER_VARIABLE( ARC_LENGTH_LAMBDA )
    KRATOS_REGISTER_VARIABLE( ARC_LENGTH_RADIUS_FACTOR )

    KRATOS_REGISTER_VARIABLE( TIME_UNIT_CONVERTER )

    KRATOS_REGISTER_VARIABLE( LOCAL_EQUIVALENT_STRAIN )
    KRATOS_REGISTER_VARIABLE( NONLOCAL_EQUIVALENT_STRAIN )

    KRATOS_REGISTER_VARIABLE( JOINT_WIDTH )

    KRATOS_REGISTER_VARIABLE( NODAL_SMOOTHING )
    KRATOS_REGISTER_VARIABLE( NODAL_CAUCHY_STRESS_TENSOR )
    KRATOS_REGISTER_VARIABLE( NODAL_DAMAGE_VARIABLE )
    KRATOS_REGISTER_VARIABLE( NODAL_JOINT_AREA )
    KRATOS_REGISTER_VARIABLE( NODAL_JOINT_WIDTH )
    KRATOS_REGISTER_VARIABLE( NODAL_JOINT_DAMAGE )

    KRATOS_REGISTER_VARIABLE( RESET_DISPLACEMENTS )
    KRATOS_REGISTER_VARIABLE( CONSIDER_GEOMETRIC_STIFFNESS )

    KRATOS_REGISTER_VARIABLE( USE_CONSISTENT_MASS_MATRIX)

    KRATOS_REGISTER_VARIABLE( IGNORE_UNDRAINED )

    KRATOS_REGISTER_VARIABLE( SATURATED_SATURATION )
    KRATOS_REGISTER_VARIABLE( RESIDUAL_SATURATION )
    KRATOS_REGISTER_VARIABLE( VAN_GENUCHTEN_AIR_ENTRY_PRESSURE )
    KRATOS_REGISTER_VARIABLE( VAN_GENUCHTEN_GN )
    KRATOS_REGISTER_VARIABLE( VAN_GENUCHTEN_GL )
    KRATOS_REGISTER_VARIABLE( MINIMUM_RELATIVE_PERMEABILITY )

    KRATOS_REGISTER_VARIABLE( RETENTION_LAW )
    KRATOS_REGISTER_VARIABLE( DEGREE_OF_SATURATION )
    KRATOS_REGISTER_VARIABLE( EFFECTIVE_SATURATION )
    KRATOS_REGISTER_VARIABLE( BISHOP_COEFICIENT )
    KRATOS_REGISTER_VARIABLE( DERIVATIVE_OF_SATURATION )
    KRATOS_REGISTER_VARIABLE( RELATIVE_PERMEABILITY )

    // UDSM
    KRATOS_REGISTER_VARIABLE( UDSM_NAME )       // Also for UMAT
    KRATOS_REGISTER_VARIABLE( UDSM_NUMBER )
    KRATOS_REGISTER_VARIABLE( IS_FORTRAN_UDSM ) // Also for UMAT

    KRATOS_REGISTER_VARIABLE( UMAT_PARAMETERS )

    KRATOS_REGISTER_VARIABLE(NUMBER_OF_UMAT_STATE_VARIABLES)
    KRATOS_REGISTER_VARIABLE(NUMBER_OF_UMAT_PARAMETERS)

    KRATOS_REGISTER_VARIABLE( STATE_VARIABLES )

    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_1 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_2 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_3 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_4 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_5 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_6 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_7 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_8 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_9 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_10 )

    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_11 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_12 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_13 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_14 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_15 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_16 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_17 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_18 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_19 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_20 )

    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_21 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_22 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_23 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_24 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_25 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_26 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_27 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_28 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_29 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_30 )

    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_31 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_32 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_33 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_34 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_35 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_36 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_37 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_38 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_39 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_40 )

    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_41 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_42 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_43 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_44 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_45 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_46 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_47 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_48 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_49 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_50 )

   }
}  // namespace Kratos.
