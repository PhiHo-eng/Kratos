//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main author:     Miguel Angel Celigueta
//

#if !defined(KRATOS_APPLY_VELOCITY_CONSTRAINTS_PROCESS )
#define  KRATOS_APPLY_VELOCITY_CONSTRAINTS_PROCESS

// System includes
#include <string>
#include <iostream>

// Project includes
#include "includes/table.h"
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "utilities/interval_utility.h"
#include "utilities/python_function_callback_utility.h"
#include "DEM_application_variables.h"

namespace Kratos
{
  ///@addtogroup DEMApplication
  ///@{

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

  /// Auxiliary process to apply forces and moments to particles.
  /** This process sets the EXTERNAL_APPLIED_FORCE and EXTERNAL_APPLIED_MOMENT variables
      over particles.
   */


  class KRATOS_API(DEM_APPLICATION) ApplyVelocityConstraintsProcess: public Process
  {
  public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of ApplyVelocityConstraintsProcess
      KRATOS_CLASS_POINTER_DEFINITION(ApplyVelocityConstraintsProcess);

      /// Defining a table with double argument and result type as table type.
      typedef Table<double,double> TableType;

      ///@}
      ///@name Life Cycle
      ///@{

      /// Constructor
      ApplyVelocityConstraintsProcess(ModelPart& rModelPart, Parameters rParameters);

      /// Destructor.
      ~ApplyVelocityConstraintsProcess() override;


      ///@}
      ///@name Operators
      ///@{


      ///@}
      ///@name Operations
      ///@{

      void Execute() override;

      void ExecuteInitialize() override;

      void ExecuteInitializeSolutionStep() override;

      void ExecuteFinalizeSolutionStep() override;

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
      std::string Info() const override;

      /// Print information about this object.
      void PrintInfo(std::ostream& rOStream) const override;

      /// Print object's data.
      void PrintData(std::ostream& rOStream) const override;


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


      ///@}
      ///@name Member Variables
      ///@{

      ModelPart& mrModelPart;
      Parameters mParameters;
      IntervalUtility mInterval;
      array_1d<bool, 3> mVelocityIsConstrained;
      array_1d<bool, 3> mVelocityValueIsNumeric;
      array_1d<double, 3> mVelocityValues;
      std::vector<GenericFunctionUtility> mVelocityFunctions;
      array_1d<int, 3> mVelocityTableId;
      std::vector<TableType::Pointer> mpVelocityTable;

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
      ApplyVelocityConstraintsProcess& operator=(ApplyVelocityConstraintsProcess const& rOther);

      /// Copy constructor.
      ApplyVelocityConstraintsProcess(ApplyVelocityConstraintsProcess const& rOther);


      ///@}

    }; // Class ApplyVelocityConstraintsProcess

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// output stream function
  inline std::ostream& operator << (
      std::ostream& rOStream,
      const ApplyVelocityConstraintsProcess& rThis);

  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_APPLY_VELOCITY_CONSTRAINTS_PROCESS  defined
