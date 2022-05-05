//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Hofer
//					 Erich Wehrle
//
// ==============================================================================

#if !defined(KRATOS_RESIDUALBASED_ADJOINTUPDATE_STATIC_SIMP_SCHEME_H)
#define  KRATOS_RESIDUALBASED_ADJOINTUPDATE_STATIC_SIMP_SCHEME_H


// External includes


// Project includes
#include "includes/define.h"

#include "topology_optimization_application.h"
#include "structural_mechanics_application.h"

// Application includes
#include "solving_strategies/schemes/residual_based_adjoint_static_scheme.h"
#include "includes/ublas_interface.h"
#include "utilities/openmp_utils.h"
#include "solving_strategies/schemes/scheme.h"
#include "response_functions/adjoint_response_function.h"


namespace Kratos
{

/**@name Kratos Globals */
/*@{ */


/*@} */
/**@name Type Definitions */
/*@{ */

/*@} */


/**@name  Enum's */
/*@{ */


/*@} */
/**@name  Functions */
/*@{ */


/*@} */
/**@name Kratos Classes */
/*@{ */

template<class TSparseSpace,class TDenseSpace >
class ResidualBasedAdjointUpdateStaticSIMPScheme : public ResidualBasedAdjointStaticScheme<TSparseSpace,TDenseSpace>
{

public:
    /**@name Type Definitions */
    /*@{ */

    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedAdjointUpdateStaticSIMPScheme);

    typedef ResidualBasedAdjointStaticScheme<TSparseSpace, TDenseSpace> BaseType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    explicit ResidualBasedAdjointUpdateStaticSIMPScheme(AdjointResponseFunction::Pointer pResponseFunction)
        : ResidualBasedAdjointStaticScheme<TSparseSpace, TDenseSpace>(pResponseFunction)
    {
        mpResponseFunction = pResponseFunction;

        int num_threads = ParallelUtilities::GetNumThreads();
        mAdjointValues.resize(num_threads);
    }

    /// Destructor.
    ~ResidualBasedAdjointUpdateStaticSIMPScheme() override
    {
    }

    /*@} */
    /**@name Operators
    */
    /*@{ */

    // ---------------------------------------------------------------------------------------------------------------------------------------------
    // --------------------------------- UPDATE LHS AND RHS ----------------------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------------------------------------------------------------------

    /// This function calculates a new Youngs Modulus based on the densities and multiplies it into the
    /// LHS and RHS contributions of the complete system



    void CalculateSystemContributions(Element& rCurrentElement,
                                      LocalSystemMatrixType& rLHS_Contribution,
                                      LocalSystemVectorType& rRHS_Contribution,
                                      Element::EquationIdVectorType& rEquationId,
                                      const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        int thread_id = OpenMPUtils::ThisThread();

        const auto& r_const_elem_ref = rCurrentElement;

        //Determine the new Youngs Modulus based on the assigned new density (X_PHYS)
        const double E_min     = rCurrentElement.GetProperties()[YOUNGS_MODULUS_MIN];
        const double E_initial = rCurrentElement.GetProperties()[YOUNGS_MODULUS_0];
        const double E_current = rCurrentElement.GetValue(YOUNG_MODULUS);
        const double penalty   = rCurrentElement.GetValue(PENAL);
        const double x_phys    = rCurrentElement.GetValue(X_PHYS);

        if(E_min == 0)
        {
            KRATOS_INFO("[TopOpt]") << "In Residual Scheme E Min: "<<E_min<< std::endl;
        }

        if(E_initial == 0)
        {
            KRATOS_INFO("[TopOpt]") << "In Residual Scheme E Initial: "<<E_initial<< std::endl;
        }

        if(E_current == 0)
        {
            KRATOS_INFO("[TopOpt]") << "In Residual Scheme E_current: "<<E_current<< std::endl;
        }

        if(penalty == 0)
        {
            KRATOS_INFO("[TopOpt]") << "In Residual Scheme penalty: "<<penalty<< std::endl;
        }


        const double E_new     = (E_min + pow(x_phys, penalty) * (E_initial - E_min));

        //Calculate the factor that needs to be multiplied on the RHS and LHS
        const double factor    = E_new/E_current;

        // Factorize LHS and RHS according SIMP approach
        // Note that when this function is called, all the contributions from the force conditions are missing.
        // I.e. RHS = -K*u_init. Hence we can directly factorize LHS and RHS to obtained the modified stiffnesses
        rLHS_Contribution *= factor;
        rRHS_Contribution *= factor;

        

        rCurrentElement.CalculateLeftHandSide(rLHS_Contribution, rCurrentProcessInfo);

        if (rRHS_Contribution.size() != rLHS_Contribution.size1())
            rRHS_Contribution.resize(rLHS_Contribution.size1(), false);

        mpResponseFunction->CalculateGradient(
            rCurrentElement, rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);

        noalias(rRHS_Contribution) = -rRHS_Contribution;

        // Calculate system contributions in residual form.
        r_const_elem_ref.GetValuesVector(mAdjointValues[thread_id]);
        noalias(rRHS_Contribution) -= prod(rLHS_Contribution, mAdjointValues[thread_id]);

        r_const_elem_ref.EquationIdVector(rEquationId, rCurrentProcessInfo);

        KRATOS_CATCH("");
    }

    void CalculateLHSContribution(Element& rCurrentElement,
                                  LocalSystemMatrixType& rLHS_Contribution,
                                  Element::EquationIdVectorType& rEquationId,
                                  const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        rCurrentElement.CalculateLeftHandSide(rLHS_Contribution, rCurrentProcessInfo);
        rCurrentElement.EquationIdVector(rEquationId, rCurrentProcessInfo);

        KRATOS_CATCH("");
    }

    void CalculateSystemContributions(Condition& rCurrentCondition,
                                      LocalSystemMatrixType& rLHS_Contribution,
                                      LocalSystemVectorType& rRHS_Contribution,
                                      Condition::EquationIdVectorType& rEquationId,
                                      const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        int thread_id = OpenMPUtils::ThisThread();
        const auto& r_const_cond_ref = rCurrentCondition;
        rCurrentCondition.CalculateLeftHandSide(rLHS_Contribution, rCurrentProcessInfo);

        if (rRHS_Contribution.size() != rLHS_Contribution.size1())
            rRHS_Contribution.resize(rLHS_Contribution.size1(), false);

        mpResponseFunction->CalculateGradient(
            rCurrentCondition, rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);

        noalias(rRHS_Contribution) = -rRHS_Contribution;

        // Calculate system contributions in residual form.
        r_const_cond_ref.GetValuesVector(mAdjointValues[thread_id]);
        noalias(rRHS_Contribution) -= prod(rLHS_Contribution, mAdjointValues[thread_id]);

        r_const_cond_ref.EquationIdVector(rEquationId, rCurrentProcessInfo);

        KRATOS_CATCH("");
    }

    void CalculateLHSContribution(Condition& rCurrentCondition,
                                  LocalSystemMatrixType& rLHS_Contribution,
                                  Condition::EquationIdVectorType& rEquationId,
                                  const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        rCurrentCondition.CalculateLeftHandSide(rLHS_Contribution, rCurrentProcessInfo);
        rCurrentCondition.EquationIdVector(rEquationId, rCurrentProcessInfo);

        KRATOS_CATCH("");
    }

  
    /*@} */
    /**@name Operations */
    /*@{ */


    /*@} */
    /**@name Access */
    /*@{ */


    /*@} */
    /**@name Inquiry */
    /*@{ */


    /*@} */
    /**@name Friends */
    /*@{ */


    /*@} */

protected:
    /**@name Protected static Member Variables */
    /*@{ */


    /*@} */
    /**@name Protected member Variables */
    /*@{ */

    AdjointResponseFunction::Pointer mpResponseFunction;
    std::vector<LocalSystemVectorType> mAdjointValues;

    /*@} */
    /**@name Protected Operators*/
    /*@{ */

    /*@} */
    /**@name Protected Operations*/
    /*@{ */


    /*@} */
    /**@name Protected  Access */
    /*@{ */


    /*@} */
    /**@name Protected Inquiry */
    /*@{ */


    /*@} */
    /**@name Protected LifeCycle */
    /*@{ */



    /*@} */

private:
    /**@name Static Member Variables */
    /*@{ */


    /*@} */
    /**@name Member Variables */
    /*@{ */

    /*@} */
    /**@name Private Operators*/
    /*@{ */


    /*@} */
    /**@name Private Operations*/
    /*@{ */


    /*@} */
    /**@name Private  Access */
    /*@{ */


    /*@} */
    /**@name Private Inquiry */
    /*@{ */


    /*@} */
    /**@name Un accessible methods */
    /*@{ */


    /*@} */

}; /* Class Scheme */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_RESIDUALBASED_ADJOINTUPDATE_STATIC_SIMP_SCHEME_H  defined */
