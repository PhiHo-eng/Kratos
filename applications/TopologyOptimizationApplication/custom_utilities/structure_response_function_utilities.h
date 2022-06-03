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

#if !defined(KRATOS_STRUCTURE_RESPONSE_FUNCTION_UTILITIES_H_INCLUDED)
#define  KRATOS_STRUCTURE_RESPONSE_FUNCTION_UTILITIES_H_INCLUDED

// System includes

// External includes

// Project includes

// Application includes
#include "topology_optimization_application.h"
#include "utilities/builtin_timer.h"


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

/// Solution utility to compute structural analysis responses.
/** Detail class definition.

 */

class StructureResponseFunctionUtilities
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of StructureResponseFunctionUtilities
    KRATOS_CLASS_POINTER_DEFINITION(StructureResponseFunctionUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    StructureResponseFunctionUtilities( ModelPart& model_part )
    : mr_structure_model_part(model_part)
    {
    }

    /// Destructor.
    virtual ~StructureResponseFunctionUtilities()
    {
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // ---------------------------------------------------------------------------------------------------------------------------------------------
    // --------------------------------- COMPUTE STRAIN ENERGY -------------------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------------------------------------------------------------------

    /// Computes the strain energy as the objective function of the optimization problem.
    double ComputeStrainEnergy()
    {
        KRATOS_TRY;

        BuiltinTimer timer;
        KRATOS_INFO("[TopOpt]") << "  Start calculating strain energy."<<std::endl;

        double Out = 0.0;
        double Global_Strain_Energy = 0.0;
        // Loop over all elements to calculate their local objective function and sum it into the global objective function (Global Strain Energy)
        for( ModelPart::ElementIterator element_i = mr_structure_model_part.ElementsBegin(); element_i!= mr_structure_model_part.ElementsEnd();
                element_i++ )
        {

            element_i->Calculate(LOCAL_STRAIN_ENERGY_COMPLIANT, Out, mr_structure_model_part.GetProcessInfo());

            Global_Strain_Energy += element_i->GetValue(LOCAL_STRAIN_ENERGY_COMPLIANT);

            //if ( std::isnan(Global_Strain_Energy))
            //{
              //  KRATOS_INFO("[TopOpt]") << "  OOOOOO Shit in structure"<<Global_Strain_Energy << std::endl;
            //}
            
        }

        KRATOS_INFO("[TopOpt]") <<  "  Strain energy calculated                  [ spent time =  " << timer.ElapsedSeconds() << " ] " << std::endl;

        // Return this obtained Global Strain Energy value as the objective function of the complete system
        return Global_Strain_Energy;

        KRATOS_CATCH("");
    }

    double ComputeDisplacement()
    {
        KRATOS_TRY;

        BuiltinTimer timer;
        KRATOS_INFO("[TopOpt]") <<"  Start the calculation of the dispalcements at the defined nodes."<<std::endl;

        array_1d<double,3> vector_L = ZeroVector(3);
        array_1d<double,3> displacement = ZeroVector(3);
        int counter = 1;
        double objective_function = 0;
        std::vector<double> output_von_mises(1);

        for( ModelPart::ElementIterator element_i = mr_structure_model_part.ElementsBegin(); element_i!= mr_structure_model_part.ElementsEnd();
                element_i++ )
        {

        element_i->CalculateOnIntegrationPoints(VON_MISES_STRESS, output_von_mises, mr_structure_model_part.GetProcessInfo());
            
        }



        for(ModelPart::NodeIterator it_node = mr_structure_model_part.NodesBegin(); it_node != mr_structure_model_part.NodesEnd(); ++it_node )
        {
            if (counter ==1) //(counter == 11381 || counter==11383 ||counter ==11393|| counter==11407 ||counter == 11421) //(counter == 177 || counter == 178 || counter == 182) //(counter == 2246 || counter == 2248)//(counter == 29382 || counter==29385 ||counter == 29394|| counter==29404 ||counter == 29416 || counter==29430  ) //
            {
                
                displacement = it_node->FastGetSolutionStepValue(DISPLACEMENT);
                objective_function +=   displacement[0];
                KRATOS_INFO("[TopOpt]") <<"  Displacement: "<<objective_function <<std::endl;
            }
            counter ++;
        }

        KRATOS_INFO("[TopOpt]") <<"  Displacement calculated                [ spent time =  " << timer.ElapsedSeconds() << " ] " << std::endl;
        return objective_function;

        KRATOS_CATCH("");
    }

    double ComputeVolumeFraction()
    {
        KRATOS_TRY;

        BuiltinTimer timer;
        KRATOS_INFO("[TopOpt]") <<"  Start calculating volume fraction."<<std::endl;

        double Global_Volume_Fraction = 0.0;
        double elemental_volume = 0.0;
        double design_variable = 0.0;
        double Total_volume = 0.0;


        // Loop over all elements to obtain their X_PHYS and know how many elements the model has
        for( ModelPart::ElementIterator element_i = mr_structure_model_part.ElementsBegin(); element_i!= mr_structure_model_part.ElementsEnd();
                element_i++ )
        {

            elemental_volume = element_i->GetValue(INITIAL_ELEMENT_SIZE);
            design_variable = element_i->GetValue(X_PHYS);
            Global_Volume_Fraction += (elemental_volume*design_variable); //
            Total_volume += elemental_volume;
        }

        // Calculate and return the Global Volume Fraction by knowing how many elements the model has
        Global_Volume_Fraction = Global_Volume_Fraction/Total_volume;
        KRATOS_INFO("[TopOpt]") <<"  Global Volume Fraction: " << Global_Volume_Fraction << std::endl;
        KRATOS_INFO("[TopOpt]") <<"  Volume fraction calculated                [ spent time =  " << timer.ElapsedSeconds() << " ] " << std::endl;
        return Global_Volume_Fraction;

        KRATOS_CATCH("");
    }

    double ComputeDisplacementControlledObjective()
    {
        KRATOS_TRY;

        BuiltinTimer timer;
        KRATOS_INFO("[TopOpt]") <<"  Start the calculation of the controlled dispalcements at the defined nodes."<<std::endl;

        array_1d<double,3> vector_L = ZeroVector(3);
        array_1d<double,3> displacement = ZeroVector(3);
        int counter = 1;
        double displacement_error_1 = 0;
        double displacement_error_2 = 0;
        double displacement_error_3 = 0;
        double displacement_error_4 = 0;
        double displacement_error_5 = 0;
        double displacement_error_6 = 0;

        std::vector<double> output_von_mises(1);

        for( ModelPart::ElementIterator element_i = mr_structure_model_part.ElementsBegin(); element_i!= mr_structure_model_part.ElementsEnd();
                element_i++ )
        {

        element_i->CalculateOnIntegrationPoints(VON_MISES_STRESS, output_von_mises, mr_structure_model_part.GetProcessInfo());
            
        }



        for(ModelPart::NodeIterator it_node = mr_structure_model_part.NodesBegin(); it_node != mr_structure_model_part.NodesEnd(); ++it_node )
        {
            //NODE number 1
            if (counter ==33715) //(counter == 11381 || counter==11383 ||counter ==11393|| counter==11407 ||counter == 11421) //(counter == 177 || counter == 178 || counter == 182) //(counter == 2246 || counter == 2248)//(counter == 29382 || counter==29385 ||counter == 29394|| counter==29404 ||counter == 29416 || counter==29430  ) //
            {
                

                displacement = it_node->FastGetSolutionStepValue(DISPLACEMENT);
                displacement_error_1 = std::sqrt(std::pow(-0.812049128806962+displacement[0],2)+std::pow(3.50266712604163+displacement[1],2)+std::pow(0+displacement[2],2));
                
                KRATOS_INFO("[TopOpt]") <<"  Displacement error 1: "<<displacement_error_1 <<std::endl;
            }

            //NODE number 2
            if (counter ==18190) //(counter == 11381 || counter==11383 ||counter ==11393|| counter==11407 ||counter == 11421) //(counter == 177 || counter == 178 || counter == 182) //(counter == 2246 || counter == 2248)//(counter == 29382 || counter==29385 ||counter == 29394|| counter==29404 ||counter == 29416 || counter==29430  ) //
            {
                

                displacement = it_node->FastGetSolutionStepValue(DISPLACEMENT);
                displacement_error_2 = std::sqrt(std::pow(-2.14620900011183+displacement[0],2)+std::pow(7.48136870876835+displacement[1],2)+std::pow(0+displacement[2],2));
                
                KRATOS_INFO("[TopOpt]") <<"  Displacement error 2: "<<displacement_error_2 <<std::endl;
            }

            //NODE number 3
            if (counter ==3369) //(counter == 11381 || counter==11383 ||counter ==11393|| counter==11407 ||counter == 11421) //(counter == 177 || counter == 178 || counter == 182) //(counter == 2246 || counter == 2248)//(counter == 29382 || counter==29385 ||counter == 29394|| counter==29404 ||counter == 29416 || counter==29430  ) //
            {
                

                displacement = it_node->FastGetSolutionStepValue(DISPLACEMENT);
                displacement_error_3 = std::sqrt(std::pow(1.05937342445258+displacement[0],2)+std::pow(0+displacement[1],2)+std::pow(0+displacement[2],2));
                
                KRATOS_INFO("[TopOpt]") <<"  Displacement error 3: "<<displacement_error_3 <<std::endl;
            }

            //NODE number 4
            if (counter ==38508) //(counter == 11381 || counter==11383 ||counter ==11393|| counter==11407 ||counter == 11421) //(counter == 177 || counter == 178 || counter == 182) //(counter == 2246 || counter == 2248)//(counter == 29382 || counter==29385 ||counter == 29394|| counter==29404 ||counter == 29416 || counter==29430  ) //
            {
                

                displacement = it_node->FastGetSolutionStepValue(DISPLACEMENT);
                displacement_error_4 = std::sqrt(std::pow(-7.47767034460666+displacement[0],2)+std::pow(11.76331913057+displacement[1],2)+std::pow(0+displacement[2],2));
                
                KRATOS_INFO("[TopOpt]") <<"  Displacement error 4: "<<displacement_error_4 <<std::endl;
            }

            //NODE number 5
            if (counter ==29977) //(counter == 11381 || counter==11383 ||counter ==11393|| counter==11407 ||counter == 11421) //(counter == 177 || counter == 178 || counter == 182) //(counter == 2246 || counter == 2248)//(counter == 29382 || counter==29385 ||counter == 29394|| counter==29404 ||counter == 29416 || counter==29430  ) //
            {
                

                displacement = it_node->FastGetSolutionStepValue(DISPLACEMENT);
                displacement_error_5 = std::sqrt(std::pow(-7.60649010985585+displacement[0],2)+std::pow(6.38600432590056+displacement[1],2)+std::pow(0+displacement[2],2));
                
                KRATOS_INFO("[TopOpt]") <<"  Displacement error 5: "<<displacement_error_5 <<std::endl;
            }

            //NODE number 6
            if (counter ==77315) //(counter == 11381 || counter==11383 ||counter ==11393|| counter==11407 ||counter == 11421) //(counter == 177 || counter == 178 || counter == 182) //(counter == 2246 || counter == 2248)//(counter == 29382 || counter==29385 ||counter == 29394|| counter==29404 ||counter == 29416 || counter==29430  ) //
            {

                displacement = it_node->FastGetSolutionStepValue(DISPLACEMENT);
                displacement_error_6 = std::sqrt(std::pow(-7.55812258995996+displacement[0],2)+std::pow(2.85005978676621+displacement[1],2)+std::pow(0+displacement[2],2));
                
                KRATOS_INFO("[TopOpt]") <<"  Displacement error 6: "<<displacement_error_6 <<std::endl;
            }
            counter ++;
        }

        double epsilon = 0.0;
        //epsilon = (1/6)*(std::pow(displacement_error_1,2)+std::pow(displacement_error_2,2)+std::pow(displacement_error_3,2)+std::pow(displacement_error_4,2)+std::pow(displacement_error_5,2)+std::pow(displacement_error_6,2));
        epsilon += std::pow(displacement_error_1,2);
        epsilon += std::pow(displacement_error_2,2);
        epsilon += std::pow(displacement_error_3,2);
        epsilon += std::pow(displacement_error_4,2);
        epsilon += std::pow(displacement_error_5,2);
        epsilon += std::pow(displacement_error_6,2);

        epsilon = (epsilon/1.0);


        KRATOS_INFO("[TopOpt]") <<"  Displacement calculated                [ spent time =  " << timer.ElapsedSeconds() << " ] " << std::endl;
        KRATOS_INFO("[TopOpt]") <<"  Epsilon ist: "<<epsilon <<std::endl;
        return epsilon;



        KRATOS_CATCH("");
    }

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
        return "StructureResponseFunctionUtilities";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "StructureResponseFunctionUtilities";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
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


    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mr_structure_model_part;

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
    //StructureResponseFunctionUtilities& operator=(StructureResponseFunctionUtilities const& rOther);

    /// Copy constructor.
    //StructureResponseFunctionUtilities(StructureResponseFunctionUtilities const& rOther);


    ///@}

}; // Class StructureResponseFunctionUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif	/* KRATOS_STRUCTURE_RESPONSE_FUNCTION_UTILITIES_H_INCLUDED */