// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//                   Geiser Armin, https://github.com/armingeiser
//
// ==============================================================================

#if !defined(KRATOS_SHAPEOPTIMIZATION_APPLICATION_H_INCLUDED )
#define  KRATOS_SHAPEOPTIMIZATION_APPLICATION_H_INCLUDED

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------

#include <string>
#include <iostream>

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/kratos_application.h"

/* ELEMENTS */
#include "custom_elements/helmholtz_element.h"
#include "custom_elements/helmholtz_vec_element.h"

/* CONDITIONS */
#include "custom_conditions/surface_filter_condition.h"

// ==============================================================================

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
	/** Detail class definition.
	*/
	class KRATOS_API(SHAPE_OPTIMIZATION_APPLICATION) KratosShapeOptimizationApplication : public KratosApplication
	{
	public:
		///@name Type Definitions
		///@{


		/// Pointer definition of KratosShapeOptimizationApplication
		KRATOS_CLASS_POINTER_DEFINITION(KratosShapeOptimizationApplication);

		///@}
		///@name Life Cycle
		///@{

		/// Default constructor.
		KratosShapeOptimizationApplication();

		/// Destructor.
		~KratosShapeOptimizationApplication() override {}


		///@}
		///@name Operators
		///@{


		///@}
		///@name Operations
		///@{

	    void Register() override;



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
		std::string Info() const override
		{
			return "KratosShapeOptimizationApplication";
		}

		/// Print information about this object.
		void PrintInfo(std::ostream& rOStream) const override
		{
			rOStream << Info();
			PrintData(rOStream);
		}

		///// Print object's data.
       void PrintData(std::ostream& rOStream) const override
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

			const HelmholtzElement mHelmholtz2D3N;
			const HelmholtzElement mHelmholtz3D4N;
			const HelmholtzElement mHelmholtz3D8N;
			const HelmholtzElement mHelmholtz3D27N;
			const HelmholtzVecElement mHelmholtzVec2D3N;
			const HelmholtzVecElement mHelmholtzVec3D4N;
			const HelmholtzVecElement mHelmholtzVec3D8N;
			const HelmholtzVecElement mHelmholtzVec3D27N; 

			/* CONDITIONS*/

			// Surface conditions
			const SurfaceFilterCondition mSurfaceFilterCondition3D3N;
			const SurfaceFilterCondition mSurfaceFilterCondition3D4N;
			const SurfaceFilterCondition mSurfaceFilterCondition3D6N;
			const SurfaceFilterCondition mSurfaceFilterCondition3D8N;
			const SurfaceFilterCondition mSurfaceFilterCondition3D9N;

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
		KratosShapeOptimizationApplication& operator=(KratosShapeOptimizationApplication const& rOther);

		/// Copy constructor.
		KratosShapeOptimizationApplication(KratosShapeOptimizationApplication const& rOther);


		///@}

	}; // Class KratosShapeOptimizationApplication

	///@}


	///@name Type Definitions
	///@{


	///@}
	///@name Input and output
	///@{

	///@}


}  // namespace Kratos.

#endif // KRATOS_SHAPEOPTIMIZATION_APPLICATION_H_INCLUDED  defined


