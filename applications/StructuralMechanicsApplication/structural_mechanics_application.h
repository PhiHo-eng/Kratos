// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if !defined(KRATOS_STRUCTURAL_MECHANICS_APPLICATION_H_INCLUDED )
#define  KRATOS_STRUCTURAL_MECHANICS_APPLICATION_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"

#include "structural_mechanics_application_variables.h"

/* ELEMENTS */
/* Adding beam element */
#include "custom_elements/small_displacement_beam_element_3D2N.hpp"

/* Adding shells and membranes elements */
#include "custom_elements/isotropic_shell_element.hpp"
#include "custom_elements/membrane_element.hpp"
#include "custom_elements/shell_thick_element_3D4N.hpp"
#include "custom_elements/shell_thin_element_3D3N.hpp"
#include "custom_elements/total_lagrangian.h"
#include "custom_elements/kinematic_linear.h"
#include "custom_elements/nodal_concentrated_element.hpp"
#include "custom_elements/SprismElement3D6N.hpp"

/* CONDITIONS */
#include "custom_conditions/point_moment_3D_condition.hpp"
#include "custom_conditions/point_torque_3D_condition.hpp"
#include "custom_conditions/base_load_condition.h"
#include "custom_conditions/line_load_condition_2d.h"
#include "custom_conditions/surface_load_condition_3d.h"

/* UTILITIES */
// Cross sections
#include "custom_utilities/shell_cross_section.hpp"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/**
 * This application features Elements, Conditions, Constitutive laws and Utilities
 * for structural analysis problems
 */
class KratosStructuralMechanicsApplication : public KratosApplication
{
public:
    ///@name Type Definitions
    ///@{


    /// Pointer definition of KratosStructuralMechanicsApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosStructuralMechanicsApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosStructuralMechanicsApplication();

    /// Destructor.
    virtual ~KratosStructuralMechanicsApplication() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    virtual void Register();



    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "KratosStructuralMechanicsApplication";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
        PrintData(rOStream);
    }

    ///// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        KRATOS_WATCH("in my application");
        KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size() );
        rOStream << "Variables:" << std::endl;
        KratosComponents<VariableData>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Elements:" << std::endl;
        KratosComponents<Element>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Conditions:" << std::endl;
        KratosComponents<Condition>().PrintData(rOStream);
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{



    //       static const ApplicationCondition  msApplicationCondition;

    ///@}
    ///@name Member Variables
    ///@{


    /* ELEMENTS */
    // Adding the beam element 
    const SmallDisplacementBeamElement3D2N   mSmallDisplacementBeamElement3D2N;

    // Adding the shells elements 
    const IsotropicShellElement  mIsotropicShellElement3D3N;
    const ShellThickElement3D4N  mShellThickElement3D4N;
    const ShellThickElement3D4N  mShellThickCorotationalElement3D4N;
    const ShellThinElement3D3N   mShellThinElement3D3N;
    const ShellThinElement3D3N   mShellThinCorotationalElement3D3N;

    // Adding the membrane element 
    const MembraneElement mMembraneElement3D3N;
    
    // Adding the SPRISM element 
    const SprismElement3D6N mSprismElement3D6N;
    
    // Adding the nodal concentrated element 
    const NodalConcentratedElement mNodalConcentratedElement2D1N;
    const NodalConcentratedElement mNodalConcentratedElement3D1N;
    
    
    /* CONDITIONS*/
    // Beam moment condition
    const PointMoment3DCondition   mPointMomentCondition3D1N;
    // Torque condition
    const PointTorque3DCondition   mPointTorqueCondition3D1N;

    //total lagrangian
    const TotalLagrangian mTotalLagrangian2D3N;
    const TotalLagrangian mTotalLagrangian2D4N;
    const TotalLagrangian mTotalLagrangian2D6N;
    const TotalLagrangian mTotalLagrangian2D8N;
    const TotalLagrangian mTotalLagrangian2D9N;
    const TotalLagrangian mTotalLagrangian3D4N;
    const TotalLagrangian mTotalLagrangian3D6N;
    const TotalLagrangian mTotalLagrangian3D8N;
    const TotalLagrangian mTotalLagrangian3D10N;
    const TotalLagrangian mTotalLagrangian3D15N;
    const TotalLagrangian mTotalLagrangian3D20N;
    const TotalLagrangian mTotalLagrangian3D27N;

    const KinematicLinear mKinematicLinear2D3N;
    const KinematicLinear mKinematicLinear2D4N;
    const KinematicLinear mKinematicLinear2D6N;
    const KinematicLinear mKinematicLinear2D8N;
    const KinematicLinear mKinematicLinear2D9N;
    const KinematicLinear mKinematicLinear3D4N;
    const KinematicLinear mKinematicLinear3D6N;
    const KinematicLinear mKinematicLinear3D8N;
    const KinematicLinear mKinematicLinear3D10N;
    const KinematicLinear mKinematicLinear3D15N;
    const KinematicLinear mKinematicLinear3D20N;
    const KinematicLinear mKinematicLinear3D27N;
    

   const LineLoadCondition2D              mLineLoadCondition2D2N;
   const LineLoadCondition2D              mLineLoadCondition2D3N;

   const SurfaceLoadCondition3D    mSurfaceLoadCondition3D3N;
   const SurfaceLoadCondition3D    mSurfaceLoadCondition3D4N;
   const SurfaceLoadCondition3D    mSurfaceLoadCondition3D6N;
   const SurfaceLoadCondition3D    mSurfaceLoadCondition3D8N;
   const SurfaceLoadCondition3D    mSurfaceLoadCondition3D9N;

    

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    KratosStructuralMechanicsApplication& operator=(KratosStructuralMechanicsApplication const& rOther);

    /// Copy constructor.
    KratosStructuralMechanicsApplication(KratosStructuralMechanicsApplication const& rOther);


    ///@}

}; // Class KratosStructuralMechanicsApplication

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_STRUCTURAL_MECHANICS_APPLICATION_H_INCLUDED  defined 


