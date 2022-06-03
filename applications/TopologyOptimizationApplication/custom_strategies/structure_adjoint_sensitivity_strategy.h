//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
//                   Octaviano Malfavón Farías
//                   Eric Gonzales
//					 Philipp Hofer
//					 Erich Wehrle
//
// ==============================================================================

#if !defined(KRATOS_STRUCTURE_ADJOINT_SENSITIVITY_STRATEGY_H_INCLUDED)
#define  KRATOS_STRUCTURE_ADJOINT_SENSITIVITY_STRATEGY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "custom_elements/small_displacement_simp_element.h"
#include "utilities/builtin_timer.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "topology_optimization_application.h"


namespace Kratos {

///@addtogroup TopologyOptimizationApplication
///@{

///@name Kratos Classes
///@{

/// Solution strategy to calculate the sensitivities.
/// Derives from the previously defined Solving Strategy

template    <class TSparseSpace, 
            class TDenseSpace, 
            class TLinearSolver
            >
class StructureAdjointSensitivityStrategy
    : public SolvingStrategy<TSparseSpace, TDenseSpace>
{
public:

    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(StructureAdjointSensitivityStrategy);

    typedef SolvingStrategy<TSparseSpace,TDenseSpace> BaseType;

    typedef typename Scheme<TSparseSpace,TDenseSpace>::Pointer SchemePointerType;
    typedef typename BuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>::Pointer BuilderAndSolverPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    StructureAdjointSensitivityStrategy( ModelPart& rStructureModelPart,
            typename TLinearSolver::Pointer pNewLinearSolver,
            const int dimension = 3)
    : BaseType(rStructureModelPart),
        mr_structure_model_part(rStructureModelPart),
        m_dimension(dimension)
    {}


    ///virtual ~StructureAdjointSensitivityStrategy()
    ~StructureAdjointSensitivityStrategy()	override
    {}

    ///@}
    ///@name Operations
    ///@{

    // ---------------------------------------------------------------------------------------------------------------------------------------------
    // --------------------------------- COMPUTE SENSITIVITIES  ------------------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------------------------------------------------------------------

    /// Computes DCDX sensitivities from the adjoint solution
    void ComputeStrainEnergySensitivities()
    {
       KRATOS_TRY;

        double Out = 0.0;
        int i= 0;

        BuiltinTimer timer;
        const ProcessInfo& ConstProcessInfo= mr_structure_model_part.GetProcessInfo();

        block_for_each(mr_structure_model_part.Elements(), [&](Element& element_i)
        {
            element_i.Calculate(DCDX, Out, ConstProcessInfo);
            i++;
            
        });

        KRATOS_INFO("[TopOpt]") << "  Objective Function sensitivities computed  [ spent time =  " << timer.ElapsedSeconds() << " ] " << std::endl;

        KRATOS_CATCH("");
    }

    void ComputeDisplacementSensitivities()
    {
       KRATOS_TRY;

        double Out = 0.0;
        int i= 0;

        BuiltinTimer timer;
        const ProcessInfo& ConstProcessInfo= mr_structure_model_part.GetProcessInfo();

        block_for_each(mr_structure_model_part.Elements(), [&](Element& element_i)
        {
            element_i.Calculate(DCDX_COMPLIANT, Out, ConstProcessInfo);
            i++;
            
        });

        KRATOS_INFO("[TopOpt]") << "  Objective Function sensitivities computed  [ spent time =  " << timer.ElapsedSeconds() << " ] " << std::endl;

        KRATOS_CATCH("");
    }


    /// Computes DVDX sensitivities from the adjoint solution
    void ComputeVolumeFractionSensitivities()
    {
        KRATOS_TRY;

        double Out = 0.0;

        BuiltinTimer timer;
        const ProcessInfo& ConstProcessInfo= mr_structure_model_part.GetProcessInfo();

        block_for_each(mr_structure_model_part.Elements(), [&](Element& element_i)
        {
            element_i.Calculate(DVDX, Out, ConstProcessInfo);
        });

        KRATOS_INFO("[TopOpt]") << "  Volume fraction sensitivities computed     [ spent time =  " << timer.ElapsedSeconds() << " ] " << std::endl;

        KRATOS_CATCH("");
    }

    void ComputeDisplacementControlledSensitivities()
    {
       KRATOS_TRY;

        double Out_1 = 0.0;
        double Out_2 = 0.0;
        double Out_3 = 0.0;

        int i= 0;
        const ProcessInfo& ConstProcessInfo= mr_structure_model_part.GetProcessInfo();

        BuiltinTimer timer;

        double displacement_error_1 = 0;
        double displacement_error_2 = 0;
        double displacement_error_3 = 0;
        double displacement_error_4 = 0;
        double displacement_error_5 = 0;
        double displacement_error_6 = 0;
        array_1d<double,3> difference_1= ZeroVector(3);
        array_1d<double,3> difference_2= ZeroVector(3);
        array_1d<double,3> difference_3= ZeroVector(3);
        array_1d<double,3> difference_4= ZeroVector(3);
        array_1d<double,3> difference_5= ZeroVector(3);
        array_1d<double,3> difference_6= ZeroVector(3);

        array_1d<double,3> displacement = ZeroVector(3);

        int counter = 0;

        for(ModelPart::NodeIterator it_node = mr_structure_model_part.NodesBegin(); it_node != mr_structure_model_part.NodesEnd(); ++it_node )
        {
            //NODE number 1
            if (counter ==33715) //(counter == 11381 || counter==11383 ||counter ==11393|| counter==11407 ||counter == 11421) //(counter == 177 || counter == 178 || counter == 182) //(counter == 2246 || counter == 2248)//(counter == 29382 || counter==29385 ||counter == 29394|| counter==29404 ||counter == 29416 || counter==29430  ) //
            {
                

                displacement = it_node->FastGetSolutionStepValue(DISPLACEMENT);
                displacement_error_1 = std::sqrt(std::pow(-0.812049128806962+displacement[0],2)+std::pow(3.50266712604163+displacement[1],2)+std::pow(0+displacement[2],2));
                difference_1[0]=-0.812049128806962+displacement[0];
                difference_1[1]=3.50266712604163+displacement[1];
                difference_1[2]=0+displacement[2];

                KRATOS_INFO("[TopOpt]") <<"  Displacement error 1: "<<displacement_error_1 <<std::endl;
            }

            //NODE number 2
            if (counter ==18190) //(counter == 11381 || counter==11383 ||counter ==11393|| counter==11407 ||counter == 11421) //(counter == 177 || counter == 178 || counter == 182) //(counter == 2246 || counter == 2248)//(counter == 29382 || counter==29385 ||counter == 29394|| counter==29404 ||counter == 29416 || counter==29430  ) //
            {
                

                displacement = it_node->FastGetSolutionStepValue(DISPLACEMENT);
                displacement_error_2 = std::sqrt(std::pow(-2.14620900011183+displacement[0],2)+std::pow(7.48136870876835+displacement[1],2)+std::pow(0+displacement[2],2));
                difference_2[0]=-2.14620900011183+displacement[0];
                difference_2[1]=7.48136870876835+displacement[1];
                difference_2[2]=0+displacement[2];

                KRATOS_INFO("[TopOpt]") <<"  Displacement error 2: "<<displacement_error_2 <<std::endl;
            }

            //NODE number 3
            if (counter ==3369) //(counter == 11381 || counter==11383 ||counter ==11393|| counter==11407 ||counter == 11421) //(counter == 177 || counter == 178 || counter == 182) //(counter == 2246 || counter == 2248)//(counter == 29382 || counter==29385 ||counter == 29394|| counter==29404 ||counter == 29416 || counter==29430  ) //
            {
                

                displacement = it_node->FastGetSolutionStepValue(DISPLACEMENT);
                displacement_error_3 = std::sqrt(std::pow(1.05937342445258+displacement[0],2)+std::pow(0+displacement[1],2)+std::pow(0+displacement[2],2));
                difference_3[0]=1.05937342445258+displacement[0];
                difference_3[1]=11.8582156396755+displacement[1];
                difference_3[2]=0+displacement[2];

                KRATOS_INFO("[TopOpt]") <<"  Displacement error 3: "<<displacement_error_3 <<std::endl;
            }

            //NODE number 4
            if (counter ==38508) //(counter == 11381 || counter==11383 ||counter ==11393|| counter==11407 ||counter == 11421) //(counter == 177 || counter == 178 || counter == 182) //(counter == 2246 || counter == 2248)//(counter == 29382 || counter==29385 ||counter == 29394|| counter==29404 ||counter == 29416 || counter==29430  ) //
            {
                

                displacement = it_node->FastGetSolutionStepValue(DISPLACEMENT);
                displacement_error_4 = std::sqrt(std::pow(-7.47767034460666+displacement[0],2)+std::pow(11.76331913057+displacement[1],2)+std::pow(0+displacement[2],2));
                difference_4[0]=-7.47767034460666+displacement[0];
                difference_4[1]=11.76331913057+displacement[1];
                difference_4[2]=0+displacement[2];

                KRATOS_INFO("[TopOpt]") <<"  Displacement error 4: "<<displacement_error_4 <<std::endl;
            }

            //NODE number 5
            if (counter ==29977) //(counter == 11381 || counter==11383 ||counter ==11393|| counter==11407 ||counter == 11421) //(counter == 177 || counter == 178 || counter == 182) //(counter == 2246 || counter == 2248)//(counter == 29382 || counter==29385 ||counter == 29394|| counter==29404 ||counter == 29416 || counter==29430  ) //
            {
                

                displacement = it_node->FastGetSolutionStepValue(DISPLACEMENT);
                displacement_error_5 = std::sqrt(std::pow(-7.60649010985585+displacement[0],2)+std::pow(6.38600432590056+displacement[1],2)+std::pow(0+displacement[2],2));
                difference_5[0]=-7.60649010985585+displacement[0];
                difference_5[1]=6.38600432590056+displacement[1];
                difference_5[2]=0+displacement[2];

                KRATOS_INFO("[TopOpt]") <<"  Displacement error 5: "<<displacement_error_5 <<std::endl;
            }

            //NODE number 6
            if (counter ==77315) //(counter == 11381 || counter==11383 ||counter ==11393|| counter==11407 ||counter == 11421) //(counter == 177 || counter == 178 || counter == 182) //(counter == 2246 || counter == 2248)//(counter == 29382 || counter==29385 ||counter == 29394|| counter==29404 ||counter == 29416 || counter==29430  ) //
            {

                displacement = it_node->FastGetSolutionStepValue(DISPLACEMENT);
                displacement_error_6 = std::sqrt(std::pow(-7.55812258995996+displacement[0],2)+std::pow(2.85005978676621+displacement[1],2)+std::pow(0+displacement[2],2));
                difference_6[0]=-7.55812258995996+displacement[0];
                difference_6[1]=2.85005978676621+displacement[1];
                difference_6[2]=0+displacement[2];

                KRATOS_INFO("[TopOpt]") <<"  Displacement error 6: "<<displacement_error_6 <<std::endl;
            }
            counter ++;
        }

        double Out = 0.0;

        for( ModelPart::ElementIterator element_i = mr_structure_model_part.ElementsBegin(); element_i!= mr_structure_model_part.ElementsEnd();
                element_i++ )
        {
            
            double sensitivity = 0.0;

            element_i->Calculate(DCDX_COMPLIANT, Out, ConstProcessInfo);
            element_i->Calculate(DCDX_COMPLIANT, Out_2, ConstProcessInfo);
            element_i->Calculate(DCDX_COMPLIANT, Out_3, ConstProcessInfo);

            double dcdx_controlled_x = element_i->GetValue(DCDX_COMPLIANT);
            double dcdx_controlled_y = element_i->GetValue(DCDX_COMPLIANT);
            double dcdx_controlled_z = element_i->GetValue(DCDX_COMPLIANT);

            sensitivity += 2*(difference_1[0]*(dcdx_controlled_x)+difference_1[1]*(dcdx_controlled_y)+difference_1[2]*(dcdx_controlled_z));
            sensitivity += 2*(difference_2[0]*(dcdx_controlled_x)+difference_2[1]*(dcdx_controlled_y)+difference_2[2]*(dcdx_controlled_z));
            sensitivity += 2*(difference_3[0]*(dcdx_controlled_x)+difference_3[1]*(dcdx_controlled_y)+difference_3[2]*(dcdx_controlled_z));
            sensitivity += 2*(difference_4[0]*(dcdx_controlled_x)+difference_4[1]*(dcdx_controlled_y)+difference_4[2]*(dcdx_controlled_z));
            sensitivity += 2*(difference_5[0]*(dcdx_controlled_x)+difference_5[1]*(dcdx_controlled_y)+difference_5[2]*(dcdx_controlled_z));
            sensitivity += 2*(difference_6[0]*(dcdx_controlled_x)+difference_6[1]*(dcdx_controlled_y)+difference_6[2]*(dcdx_controlled_z));

            sensitivity *= (2.0/1.0);
            //KRATOS_INFO("[TopOpt]") <<"  Sensitivity: "<<sensitivity <<std::endl;
            element_i->SetValue(DCDX_COMPLIANT, sensitivity);

            
        }

        KRATOS_INFO("[TopOpt]") << "  Objective Function sensitivities computed  [ spent time =  " << timer.ElapsedSeconds() << " ] " << std::endl;

        KRATOS_CATCH("");
    }

    ///@}

private:

    ///@name Member Variables
    ///@{

    ModelPart& mr_structure_model_part;
    ModelPart* mpAdjointModelPart;
    typename BaseType::Pointer mpStrategy;
    int m_dimension;

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
}; // class StructureAdjointSensitivityStrategy

///@} // Kratos classes
///@} // AdjointStructureApplication group
}

#endif	/* KRATOS_STRUCTURE_ADJOINT_SENSITIVITY_STRATEGY_H_INCLUDED */