// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ / 
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /  
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:         BSD License
//                   license: ContactStructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"
#include "includes/mapping_variables.h"
#include "includes/master_slave_constraint.h"
#include "custom_frictional_laws/frictional_law_with_derivative.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    using GeometryType = Geometry<Node>;

///@}
///@name  Enum's
///@{

    /**
     * @brief An enumeration of the different cases of friction that can be used in the contact search process.
     */
    enum class FrictionalCase {
        FRICTIONLESS = 0, ///< Frictionless contact.
        FRICTIONLESS_COMPONENTS = 1, ///< Frictionless contact with components.
        FRICTIONAL = 2, ///< Frictional contact.
        FRICTIONLESS_PENALTY = 3, ///< Frictionless contact with penalty method.
        FRICTIONAL_PENALTY = 4 ///< Frictional contact with penalty method.
    };

    /**
     * @brief An enumeration of the different options for normal derivatives computation.
     */
    enum NormalDerivativesComputation {
        NO_DERIVATIVES_COMPUTATION = 0, ///< No computation of normal derivatives.
        ELEMENTAL_DERIVATIVES = 1, ///< Computation of normal derivatives at the element level.
        NODAL_ELEMENTAL_DERIVATIVES = 2, ///< Computation of normal derivatives at the node level.
        NO_DERIVATIVES_COMPUTATION_WITH_NORMAL_UPDATE = 3 ///< No computation of normal derivatives, but update of the normal vector.
    };

///@}
///@name  Functions
///@{

// VARIABLES
// MPC Contact related variables
KRATOS_DEFINE_APPLICATION_VARIABLE( CONTACT_STRUCTURAL_MECHANICS_APPLICATION, MasterSlaveConstraint::Pointer, CONSTRAINT_POINTER )     // Pointer to the constraint of the condition
KRATOS_DEFINE_APPLICATION_VARIABLE( CONTACT_STRUCTURAL_MECHANICS_APPLICATION, double, REACTION_CHECK_STIFFNESS_FACTOR )                // The reaction factor to be considered on the tension check

/* Mortar method */
KRATOS_DEFINE_APPLICATION_VARIABLE( CONTACT_STRUCTURAL_MECHANICS_APPLICATION, bool, CONSIDER_TESSELLATION )                            // If we consider tesellation when doing the mortar segmentation
KRATOS_DEFINE_APPLICATION_VARIABLE( CONTACT_STRUCTURAL_MECHANICS_APPLICATION, int , INNER_LOOP_ITERATION )                             // The number of loops in the simplified semi-smooth inner iteration
KRATOS_DEFINE_APPLICATION_VARIABLE( CONTACT_STRUCTURAL_MECHANICS_APPLICATION, int , INTEGRATION_ORDER_CONTACT )                        // The integration order computed in the contact
KRATOS_DEFINE_APPLICATION_VARIABLE( CONTACT_STRUCTURAL_MECHANICS_APPLICATION, double, DISTANCE_THRESHOLD )                             // The distance threshold considered
KRATOS_DEFINE_APPLICATION_VARIABLE( CONTACT_STRUCTURAL_MECHANICS_APPLICATION, double, ZERO_TOLERANCE_FACTOR )                          // The epsilon factor considered
KRATOS_DEFINE_APPLICATION_VARIABLE( CONTACT_STRUCTURAL_MECHANICS_APPLICATION, double, ACTIVE_CHECK_FACTOR )                            // The factor employed to search an active/inactive node
KRATOS_DEFINE_APPLICATION_VARIABLE( CONTACT_STRUCTURAL_MECHANICS_APPLICATION, double, SLIP_THRESHOLD )                                 // The threshold employed to search an slip/stick node
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(CONTACT_STRUCTURAL_MECHANICS_APPLICATION, AUXILIAR_COORDINATES )                 // Auxiliary coordinates used to map
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(CONTACT_STRUCTURAL_MECHANICS_APPLICATION, DELTA_COORDINATES )                    // Delta coordinates used to map
KRATOS_DEFINE_APPLICATION_VARIABLE( CONTACT_STRUCTURAL_MECHANICS_APPLICATION, double, NORMAL_GAP )                                     // The normal gap employed in contact formulation
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(CONTACT_STRUCTURAL_MECHANICS_APPLICATION, TANGENT_SLIP )                         // The tangent slip gap employed in contact formulation

/* Weighted values */
KRATOS_DEFINE_APPLICATION_VARIABLE( CONTACT_STRUCTURAL_MECHANICS_APPLICATION, double, WEIGHTED_GAP )                                   // The integrated gap employed in mortar formulation
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(CONTACT_STRUCTURAL_MECHANICS_APPLICATION, WEIGHTED_SLIP )                        // The integrated slip employed in mortar formulation
KRATOS_DEFINE_APPLICATION_VARIABLE( CONTACT_STRUCTURAL_MECHANICS_APPLICATION, double, WEIGHTED_SCALAR_RESIDUAL )                       // The integrated scalar residual
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(CONTACT_STRUCTURAL_MECHANICS_APPLICATION, WEIGHTED_VECTOR_RESIDUAL )             // The integrated vector residual

/* For ALM mortar condition */
KRATOS_DEFINE_APPLICATION_VARIABLE( CONTACT_STRUCTURAL_MECHANICS_APPLICATION, bool, ACTIVE_SET_COMPUTED )                              // To know if the active set has been computed
KRATOS_DEFINE_APPLICATION_VARIABLE( CONTACT_STRUCTURAL_MECHANICS_APPLICATION, bool, ACTIVE_SET_CONVERGED )                             // To know if the active set has converged
KRATOS_DEFINE_APPLICATION_VARIABLE( CONTACT_STRUCTURAL_MECHANICS_APPLICATION, bool, SLIP_SET_CONVERGED )                               // To know if the slip set has converged
KRATOS_DEFINE_APPLICATION_VARIABLE( CONTACT_STRUCTURAL_MECHANICS_APPLICATION, double, OPERATOR_THRESHOLD )                             // Consider objetive/non-objetive formulation threshold
KRATOS_DEFINE_APPLICATION_VARIABLE( CONTACT_STRUCTURAL_MECHANICS_APPLICATION, double, SLIP_AUGMENTATION_COEFFICIENT )                  // Coefficient to improve the slip computation convergence (augmented part related)
KRATOS_DEFINE_APPLICATION_VARIABLE( CONTACT_STRUCTURAL_MECHANICS_APPLICATION, double, DYNAMIC_FACTOR )                                 // The factor considered for dynamic problems (in order to take intro account the gap evolution)
KRATOS_DEFINE_APPLICATION_VARIABLE( CONTACT_STRUCTURAL_MECHANICS_APPLICATION, double, LAGRANGE_MULTIPLIER_CONTACT_PRESSURE )           // The lagrange multiplier for normal contact pressure
KRATOS_DEFINE_APPLICATION_VARIABLE( CONTACT_STRUCTURAL_MECHANICS_APPLICATION, double, AUGMENTED_NORMAL_CONTACT_PRESSURE )              // The resultant augmented pressure in the normal direction
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(CONTACT_STRUCTURAL_MECHANICS_APPLICATION, AUGMENTED_TANGENT_CONTACT_PRESSURE )   // The resultant augmented pressure in the tangent direction
KRATOS_DEFINE_APPLICATION_VARIABLE( CONTACT_STRUCTURAL_MECHANICS_APPLICATION, int, CONSIDER_NORMAL_VARIATION )                         // A value used to check if consider normal variation or not
KRATOS_DEFINE_APPLICATION_VARIABLE( CONTACT_STRUCTURAL_MECHANICS_APPLICATION, bool, ADAPT_PENALTY )                                    // To set if the penalty is recalculated or not
KRATOS_DEFINE_APPLICATION_VARIABLE( CONTACT_STRUCTURAL_MECHANICS_APPLICATION, double, MAX_GAP_FACTOR )                                 // The factor between the nodal H and the max gap considered to recalculate the penalty

/* For mesh tying mortar condition */
KRATOS_DEFINE_APPLICATION_VARIABLE( CONTACT_STRUCTURAL_MECHANICS_APPLICATION, std::string, TYING_VARIABLE )                            // The variable name for the mesh tying

/* Explicit simulation */
KRATOS_DEFINE_APPLICATION_VARIABLE( CONTACT_STRUCTURAL_MECHANICS_APPLICATION, double, MAX_GAP_THRESHOLD )                              // The gap considered as threshold to rescale penalty

/* Frictional laws */
KRATOS_DEFINE_APPLICATION_VARIABLE( CONTACT_STRUCTURAL_MECHANICS_APPLICATION, FrictionalLaw::Pointer, FRICTIONAL_LAW )                 // The frictional law considered
KRATOS_DEFINE_APPLICATION_VARIABLE( CONTACT_STRUCTURAL_MECHANICS_APPLICATION, double, TRESCA_FRICTION_THRESHOLD )                      // The threshold value for Tresca frictional contact
}