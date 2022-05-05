// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    E. Wehrle
//                   based on max stress nad previous PR to aggregation from Baumgaertner Daniel, https://github.com/dbaumgaertner
//

// System includes

// External includes

// Project includes
#include "adjoint_KS_max_stress_response_function.h"
#include "utilities/compare_elements_and_conditions_utility.h"

namespace Kratos
{
    AdjointKSMaxStressResponseFunction::AdjointKSMaxStressResponseFunction(ModelPart& rAdjointModelPart, Parameters ResponseSettings)
    : AdjointStructuralResponseFunction(rAdjointModelPart, ResponseSettings),
      mrAdjointModelPart(rAdjointModelPart),
      mCriticalPartName(ResponseSettings["critical_part_name"].GetString())
    {
        mTracedStressType = StressResponseDefinitions::ConvertStringToTracedStressType(ResponseSettings["stress_type"].GetString());
        mStressTreatment = StressResponseDefinitions::ConvertStringToStressTreatment(ResponseSettings["stress_treatment"].GetString());

        //pKS = ResponseSettings["aggregation_penalty"].GetDouble();

        if(ResponseSettings.Has("echo_level"))
            mEchoLevel = ResponseSettings["echo_level"].GetInt();

        KRATOS_ERROR_IF(mStressTreatment != StressTreatment::Mean)
            << "AdjointKSMaxStressResponseFunction::AdjointKSMaxStressResponseFunction: Specified stress treatment not supported: " << ResponseSettings["stress_type"].GetString() << std::endl;

        ModelPart& adjoint_aggregation_part = mrAdjointModelPart.GetSubModelPart(mCriticalPartName);
        for(auto& elem : adjoint_aggregation_part.Elements())
            mAggregatedElementIds.push_back(elem.Id());

        // Initialize
        for(auto& elem : adjoint_aggregation_part.Elements())
        {
            elem.SetValue(TRACED_STRESS_TYPE, static_cast<int>(mTracedStressType));
        }
    }

    // primal analysis for max stress approximated via modified Kreisselmeier-Steinhauser aggregation
    AdjointKSMaxStressResponseFunction::~AdjointKSMaxStressResponseFunction(){}

    double AdjointKSMaxStressResponseFunction::CalculateValue(ModelPart& rPrimalModelPart)
    {
        KRATOS_TRY;

        ModelPart& primal_agglomeration_part = rPrimalModelPart.GetSubModelPart(mCriticalPartName);

        IndexType elem_id_at_max = 1;

        max_mean_stress = 0.0;
        double min_stress = 1000000.0;
        double x_phys = 0.0;
        double x_phys_max = 0.0;

        for(auto& elem : primal_agglomeration_part.Elements())
        {
            Vector element_stress;
            const ProcessInfo &r_process_info = rPrimalModelPart.GetProcessInfo();
            StressCalculation::CalculateStressOnGP(elem,
                                                   mTracedStressType,
                                                   element_stress,
                                                   r_process_info);

            const SizeType stress_vec_size = element_stress.size();
            double mean_stress = 0.0;

            x_phys = elem.GetValue(X_PHYS);
            
            for(IndexType i = 0; i < stress_vec_size; ++i)
                mean_stress += element_stress[i]*std::pow(x_phys,q_relaxation);
            mean_stress /= stress_vec_size;

            if(mean_stress > max_mean_stress)
            {
                max_mean_stress = mean_stress;
                elem_id_at_max = elem.Id();
                x_phys_max = x_phys;
                
            }
            mean_stress_vector[elem.Id()] = mean_stress;
            
            if (mean_stress<min_stress)
                min_stress = mean_stress;
        }

        KRATOS_INFO("[TopOpt]") << "  Max stress: "<< max_mean_stress << std::endl;
        KRATOS_INFO("[TopOpt]") << "  Min Stress "<< min_stress << std::endl;
        KRATOS_INFO("[TopOpt]") << "  X_phys für max stress: "<< x_phys_max << std::endl;

        // reloop over elements, faster as vector operation?
        for(auto& elem : primal_agglomeration_part.Elements())
        {
            KS_exp_sum += std::exp(pKS*(mean_stress_vector[elem.Id()]-max_mean_stress));
        }
        double KS_max_mean_stress = max_mean_stress+1/pKS*std::log(KS_exp_sum);

        KRATOS_INFO_IF("AdjointKSMaxStressResponseFunction::CalculateValue", mEchoLevel > 0) << "Id of element with max stress value = " << elem_id_at_max << std::endl;
        KRATOS_INFO_IF("AdjointKSMaxStressResponseFunction::CalculateValue", mEchoLevel > 0) << "max mean stress approximated with KS aggregation = " << KS_max_mean_stress << std::endl;

        // Set traced element in adjoint model part corresponding to the found element in the primal part
        mpTracedElementInAdjointPart = mrAdjointModelPart.pGetElement(elem_id_at_max);
        mpTracedElementInAdjointPart->SetValue(TRACED_STRESS_TYPE, static_cast<int>(mTracedStressType) );

        return KS_max_mean_stress;

        KRATOS_CATCH("");
    }

    // adjoint equation for max stress approximated via modified Kreisselmeier-Steinhauser aggregation
    void AdjointKSMaxStressResponseFunction::CalculateGradient(const Element& rAdjointElement,
                                                               const Matrix& rResidualGradient,
                                                               Vector& rResponseGradient,
                                                               const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        //KRATOS_ERROR_IF(mpTracedElementInAdjointPart == nullptr)
        //    << "AdjointKSMaxStressResponseFunction::CalculateGradient: Traced element not initialized. First call \"CalculateValue\"!" << std::endl;


        rResponseGradient.clear();

        if(std::find(mAggregatedElementIds.begin(), mAggregatedElementIds.end(), rAdjointElement.Id()) != mAggregatedElementIds.end())
        {

            Matrix stress_displacement_derivative;

            Element::Pointer pElem = mrAdjointModelPart.pGetElement(rAdjointElement.Id());
            pElem->Calculate(STRESS_DISP_DERIV_ON_GP,
                             stress_displacement_derivative,
                             rProcessInfo);
            int iteration = 0;
            int num_of_derivatives_per_stress = stress_displacement_derivative.size1();
            int num_of_stress_positions = stress_displacement_derivative.size2();
            
            for (IndexType deriv_it = 0 ; deriv_it < num_of_derivatives_per_stress; ++deriv_it)
            {
                for(IndexType stress_it = 0; stress_it < num_of_stress_positions; ++stress_it)
                {
                    
                    //Check if nan
                    if ( std::isnan(stress_displacement_derivative(deriv_it, stress_it)))
                    {
                        iteration = 1;
                        KRATOS_INFO("[TopOpt]") << "pks parameter is (in step I): "<<pKS << std::endl;
                        KRATOS_INFO("[TopOpt]") << "max mean stress is (in step I): "<<max_mean_stress << std::endl;
                        KRATOS_INFO("[TopOpt]") << "KS approximated stress is (in step I): "<<KS_exp_sum << std::endl;
                        KRATOS_INFO("[TopOpt]") << "q_relaxiation is (in step I): "<<q_relaxation << std::endl;
                        KRATOS_INFO("[TopOpt]") << "Number of iteration (in step I): "<<iteration<< std::endl;
                        if (rAdjointElement.Id()==17850)
                        {
                            KRATOS_ERROR << "The value of the gradient is nan and the element Id is: " <<stress_displacement_derivative <<std::endl;
                        }

                        KRATOS_INFO("[TopOpt]") << "The value of the gradient is nan and the element Id is: " <<deriv_it<<" and "<< stress_it <<std::endl;
                    }
                }
            }
            if (iteration==1)
            {
                KRATOS_INFO("[TopOpt]") << "STRESS SENSITIVITY DERIVATIVE (in step I): "<<stress_displacement_derivative<< std::endl;
            }
            this->ExtractMeanStressDerivative(stress_displacement_derivative,
                                              rResponseGradient);

            KRATOS_ERROR_IF(rResponseGradient.size() != rResidualGradient.size1())
                << "AdjointKSMaxStressResponseFunction::CalculateGradient: Size of stress displacement derivative does not fit!" << std::endl;

            //double mean_stress = mean_stress_vector[rAdjointElement.Id()];
            double x_phys = rAdjointElement.GetValue(X_PHYS);
            //KRATOS_INFO("[TopOpt]") << "  \nResponse Gradient rResponseGradient with ID: "<<rAdjointElement.Id() <<" is:\n "<<rResponseGradient<< std::endl;
            //KRATOS_INFO("[TopOpt]") << "  \nThe stress displacement derivative 'stress_displacement_derivative' with ID: "<<rAdjointElement.Id() <<" is:\n "<<stress_displacement_derivative << std::endl;
            rResponseGradient *= ((-1) *std::exp(pKS*(mean_stress_vector[rAdjointElement.Id()]-max_mean_stress))/KS_exp_sum)*std::pow(x_phys, q_relaxation);
            int i = 0;
            int j = 0;
            for (i=0; i<rResidualGradient.size1();i++)
            {
                for (j=0; j<rResidualGradient.size2(); j++)
                {
                    if ( std::isnan(rResidualGradient(i,j)))
                    {
                        KRATOS_INFO("[TopOpt]") << "The Sensitivity Matrix is (in step I): "<<rResidualGradient << std::endl;
                        KRATOS_INFO("[TopOpt]") << "The Sensetivity Gradient is (in step I): "<<rResponseGradient << std::endl;
                        KRATOS_INFO("[TopOpt]") << "Mean stress is (in step I): "<<mean_stress_vector[rAdjointElement.Id()] << std::endl;
                        KRATOS_INFO("[TopOpt]") << "pks parameter is (in step I): "<<pKS << std::endl;
                        KRATOS_INFO("[TopOpt]") << "max mean stress is (in step I): "<<max_mean_stress << std::endl;
                        KRATOS_INFO("[TopOpt]") << "x_pys is (in step I): "<<x_phys << std::endl;
                        KRATOS_INFO("[TopOpt]") << "KS approximated stress is (in step I): "<<KS_exp_sum << std::endl;
                        KRATOS_INFO("[TopOpt]") << "q_relaxiation is (in step I): "<<q_relaxation << std::endl;
                        KRATOS_INFO("[TopOpt]") << "Column line is (in step I): "<< i <<" and "<< j << std::endl;

                        KRATOS_ERROR << "The value of the gradient is nan and the element Id is: " <<rAdjointElement.Id() <<std::endl;
                    }
                }
            }
        }
        else
        {
            rResponseGradient.clear();
        }

        KRATOS_CATCH("");
    }

    void AdjointKSMaxStressResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                                                         const Variable<double>& rVariable,
                                                                         const Matrix& rSensitivityMatrix,
                                                                         Vector& rSensitivityGradient,
                                                                         const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());

        if(std::find(mAggregatedElementIds.begin(), mAggregatedElementIds.end(), rAdjointElement.Id()) != mAggregatedElementIds.end())
        {
            this->CalculateElementContributionToPartialSensitivity(rAdjointElement,
                                                                   rVariable.Name(),
                                                                   rSensitivityMatrix,
                                                                   rSensitivityGradient,
                                                                   rProcessInfo);
            double x_phys = rAdjointElement.GetValue(X_PHYS);
             //KRATOS_INFO("[TopOpt]") << "  Gradient hier 2: " << std::endl;
            //if (rAdjointElement.Id()==100)
            //KRATOS_INFO("[TopOpt]") << "  \nThe Sensitivity Gradient rSensitivityGradient for ID : "<<rAdjointElement.Id()<<" is:\n "<< rSensitivityGradient << std::endl;
            //KRATOS_INFO("[TopOpt]") << "  \nThe Sensitivity Matrix rSensitivityMatrix for ID : "<<rAdjointElement.Id()<<" is:\n "<< rSensitivityMatrix << std::endl;
            //KRATOS_INFO("[TopOpt]") << "  \nThe Variable Name for ID : "<<rAdjointElement.Id()<<" is: "<< rVariable.Name() << std::endl;
            rSensitivityGradient *= (std::exp(pKS*(mean_stress_vector[rAdjointElement.Id()]-max_mean_stress))/KS_exp_sum)*mean_stress_vector[rAdjointElement.Id()]*q_relaxation*std::pow(x_phys, q_relaxation-1);
            //KRATOS_INFO("[TopOpt]") << "  \nThe Sensitivity Gradient rSensitivityGradient after KS for ID : "<<rAdjointElement.Id()<<" is:\n "<< rSensitivityGradient << std::endl;
            double nan_variable = rSensitivityGradient[0];
            int i = 0;
            for (i=0; i<rSensitivityGradient.size();i++)
            {
                if ( std::isnan(rSensitivityGradient[i]))
                {
                    KRATOS_INFO("[TopOpt]") << "The Sensitivity Matrix is(in step II): "<<rSensitivityMatrix << std::endl;
                    KRATOS_INFO("[TopOpt]") << "The Sensetivity Gradient is (in step II): "<<rSensitivityGradient << std::endl;
                    KRATOS_INFO("[TopOpt]") << "Mean stress is (in step II): "<<mean_stress_vector[rAdjointElement.Id()] << std::endl;
                    KRATOS_INFO("[TopOpt]") << "pks parameter is (in step II): "<<pKS << std::endl;
                    KRATOS_INFO("[TopOpt]") << "max mean stress is (in step II): "<<max_mean_stress << std::endl;
                    KRATOS_INFO("[TopOpt]") << "x_pys is (in step II): "<<x_phys << std::endl;
                    KRATOS_INFO("[TopOpt]") << "KS approximated stress is (in step II): "<<KS_exp_sum << std::endl;
                    KRATOS_INFO("[TopOpt]") << "q_relaxiation is (in step II): "<<q_relaxation << std::endl;

                    KRATOS_ERROR << "The value of the gradient is nan and the element Id is: " <<rAdjointElement.Id() <<std::endl;
                }
            }
            if ( std::isnan(nan_variable))
                {
                    KRATOS_INFO("[TopOpt]") << "  \nThe value of the gradient is nan and the element Id is : "<<rAdjointElement.Id() << std::endl;
                }
        }
        else
            rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
            //KRATOS_INFO("[TopOpt]") << "  \nThe Sensitivity Gradient rSensitivityGradient AFTER ALL SHOULD BE FINISHED for ID : "<<rAdjointElement.Id()<<" is:\n "<< rSensitivityGradient << std::endl;

        KRATOS_CATCH("")
    }

    void AdjointKSMaxStressResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                                                         const Variable<double>& rVariable,
                                                                         const Matrix& rSensitivityMatrix,
                                                                         Vector& rSensitivityGradient,
                                                                         const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());

        KRATOS_CATCH("");
    }

    void AdjointKSMaxStressResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                                                         const Variable<array_1d<double,
                                                                         3>>& rVariable,
                                                                         const Matrix& rSensitivityMatrix,
                                                                         Vector& rSensitivityGradient,
                                                                         const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;


         if(std::find(mAggregatedElementIds.begin(), mAggregatedElementIds.end(), rAdjointElement.Id()) != mAggregatedElementIds.end())
        {
            this->CalculateElementContributionToPartialSensitivity(rAdjointElement,
                                                                   rVariable.Name(),
                                                                   rSensitivityMatrix,
                                                                   rSensitivityGradient,
                                                                   rProcessInfo);

            double x_phys = rAdjointElement.GetValue(X_PHYS);
            //KRATOS_INFO("[TopOpt]") << "  Gradient hier 3: " <<rSensitivityGradient << std::endl;
            //KRATOS_INFO("[TopOpt]") << "  \nThe Sensitivity Gradient rSensitivityGradient for ID : "<<rAdjointElement.Id()<<" is:\n "<< rSensitivityGradient << std::endl;
            //KRATOS_INFO("[TopOpt]") << "  \nThe Variable name rVariable.Name() for ID : "<<rAdjointElement.Id()<<" is: "<< rVariable.Name() << std::endl;

            rSensitivityGradient *= (std::exp(pKS*(mean_stress_vector[rAdjointElement.Id()]-max_mean_stress))/KS_exp_sum)*mean_stress_vector[rAdjointElement.Id()]*q_relaxation*std::pow(x_phys,q_relaxation-1);
            int i = 0;
            for (i=0; i<rSensitivityGradient.size();i++)
            {
                if ( std::isnan(rSensitivityGradient[i]))
                {
                    KRATOS_INFO("[TopOpt]") << "The Sensitivity Matrix is (in step III): "<<rSensitivityMatrix << std::endl;
                    KRATOS_INFO("[TopOpt]") << "The Sensetivity Gradient is (in step III): "<<rSensitivityGradient << std::endl;
                    KRATOS_INFO("[TopOpt]") << "Mean stress is (in step III): "<<mean_stress_vector[rAdjointElement.Id()] << std::endl;
                    KRATOS_INFO("[TopOpt]") << "pks parameter is (in step III): "<<pKS << std::endl;
                    KRATOS_INFO("[TopOpt]") << "max mean stress is (in step III): "<<max_mean_stress << std::endl;
                    KRATOS_INFO("[TopOpt]") << "x_pys is (in step III): "<<x_phys << std::endl;
                    KRATOS_INFO("[TopOpt]") << "KS approximated stress is (in step III): "<<KS_exp_sum << std::endl;
                    KRATOS_INFO("[TopOpt]") << "q_relaxiation is (in step III): "<<q_relaxation << std::endl;

                    KRATOS_ERROR << "The value of the gradient is nan and the element Id is: " <<rAdjointElement.Id() <<std::endl;
                }
            }
        }
        else
            rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());

        KRATOS_CATCH("");
    }

    void AdjointKSMaxStressResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                                                         const Variable<array_1d<double,
                                                                         3>>& rVariable,
                                                                         const Matrix& rSensitivityMatrix,
                                                                         Vector& rSensitivityGradient,
                                                                         const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());

        KRATOS_CATCH("");
    }

    void AdjointKSMaxStressResponseFunction::CalculateElementContributionToPartialSensitivity(Element& rAdjointElement,
                                                                                              const std::string& rVariableName,
                                                                                              const Matrix& rSensitivityMatrix,
                                                                                              Vector& rSensitivityGradient,
                                                                                              const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        rAdjointElement.SetValue(DESIGN_VARIABLE_NAME, rVariableName);

        Matrix stress_design_variable_derivative;

        rAdjointElement.Calculate(STRESS_DESIGN_DERIVATIVE_ON_GP,
                                  stress_design_variable_derivative,
                                  rProcessInfo);

        this->ExtractMeanStressDerivative(stress_design_variable_derivative,
                                          rSensitivityGradient);

        KRATOS_ERROR_IF(rSensitivityGradient.size() != rSensitivityMatrix.size1()) << "Size of partial stress design variable derivative does not fit!" << std::endl;

        rAdjointElement.SetValue(DESIGN_VARIABLE_NAME, "");

        KRATOS_CATCH("");
    }

    void AdjointKSMaxStressResponseFunction::ExtractMeanStressDerivative(const Matrix& rStressDerivativesMatrix,
                                                                         Vector& rResponseGradient)
    {
        KRATOS_TRY;

        const SizeType num_of_derivatives_per_stress = rStressDerivativesMatrix.size1();
        const SizeType num_of_stress_positions = rStressDerivativesMatrix.size2();
        double stress_derivative_value = 0.0;

        if(rResponseGradient.size() != num_of_derivatives_per_stress)
            rResponseGradient.resize(num_of_derivatives_per_stress, false);

        for (IndexType deriv_it = 0 ; deriv_it < num_of_derivatives_per_stress; ++deriv_it)
        {
            for(IndexType stress_it = 0; stress_it < num_of_stress_positions; ++stress_it)
            {
                stress_derivative_value += rStressDerivativesMatrix(deriv_it, stress_it);
                
                //Check if nan
                if ( std::isnan(rStressDerivativesMatrix(deriv_it, stress_it)))
                {
                    KRATOS_INFO("[TopOpt]") << "pks parameter is (in step IIII): "<<pKS << std::endl;
                    KRATOS_INFO("[TopOpt]") << "max mean stress is (in step IIII): "<<max_mean_stress << std::endl;
                    KRATOS_INFO("[TopOpt]") << "KS approximated stress is (in step IIII): "<<KS_exp_sum << std::endl;
                    KRATOS_INFO("[TopOpt]") << "q_relaxiation is (in step IIII): "<<q_relaxation << std::endl;
                    KRATOS_INFO("[TopOpt]") << "Number of stress positions (in step IIII): "<<num_of_stress_positions << std::endl;

                    KRATOS_INFO("[TopOpt]") << "The value of the gradient is nan and the element Id is: " <<rStressDerivativesMatrix(deriv_it, stress_it) <<std::endl;
                }
            }
            stress_derivative_value /= num_of_stress_positions;

            rResponseGradient[deriv_it] = stress_derivative_value;
            stress_derivative_value = 0.0;
        }

        KRATOS_CATCH("");
    }

} // namespace Kratos.