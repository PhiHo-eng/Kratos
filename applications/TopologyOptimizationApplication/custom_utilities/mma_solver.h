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
#if !defined(KRATOS_MMA_ALGORITHM_H_INCLUDED)
#define KRATOS_MMA_ALGORITHM_H_INCLUDED

#include "includes/define.h"
#include "external_libraries/MMASolver.cpp"
#include "topology_optimization_application.h"
#include "structure_response_function_utilities.h"

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



class MMAAlgorithm 

{
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of IOUtilities
    KRATOS_CLASS_POINTER_DEFINITION(MMAAlgorithm);

	///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MMAAlgorithm(ModelPart& model_part,int n, int m)
    : mrModelPart(model_part)
    {
        this->mpMMA = new MMASolver(n,m,0,1000,0);

    }

    virtual ~MMAAlgorithm()
    {
        
    }

    void UpdateDensitiesUsingMMAAlgorithm( char update_type[], char optimization_type[], double volfrac, double greyscale , double OptItr , double qmax)
    {
        KRATOS_TRY;

        if ( strcmp( update_type , "mma_algorithm" ) == 0 )
		{
			clock_t begin = clock();
			KRATOS_INFO("[TopOpt]") << "  Method of Moving Asymthodes (MMA) chosen to solve   the optimization problem" << std::endl;


			// Create object of updating function
			int nn = mrModelPart.NumberOfElements();
			int mm = 2; /// constraints

			// number of iterations while running through the elements calling their properties
			int iteration = 0; 

			// constructing vectors for the density of the given iteration, of the iteration before and the iteration before that
			// additionally the sensitivities and constraints are saved to pass them on to the MMA optimization algorithm
            double *x = new double[nn];
			double *df = new double[nn];
			double *g = new double[mm];
			double *dg = new double[nn*mm];
			double *xmin = new double[nn];
			double *xmax = new double[nn];
			double vol_summ = 0;
			double vol_frac_iteration = 0;
            double domain_size = 0.0;
            double max_mean_stress = 0.0;
            int solid_void = 0;
            double move = 1;//0.05 changed 26.04.
            double max_stress_sens = 0.0;
            double max_vol_constraint = 0.0;
            double max_obj_func = 0.0;
            double min_stress_sens = 0.0;

			//get the information from the element and save them in the new vectors

            for( ModelPart::ElementsContainerType::iterator element_i = mrModelPart.ElementsBegin(); element_i!= mrModelPart.ElementsEnd(); element_i++ )
            {
                double initial_size = element_i->GetValue(INITIAL_ELEMENT_SIZE);

                double max_stress_sensitivity = element_i->GetValue(YOUNGS_MODULUS_SENSITIVITY);
                //KRATOS_INFO("[TopOpt]") << "  \nSENSITIVITY: "<<max_stress_sensitivity<< std::endl;

                //Value of the objective function sensitivity
                double dfdx = (element_i->GetValue(DCDX_COMPLIANT));

                // Value of the constraint function sensitivity
                double dgdx = (element_i->GetValue(DVDX));

                domain_size += initial_size;

                if (dfdx > max_obj_func)
                {
                    max_obj_func = dfdx;
                }         

                if (dgdx > max_vol_constraint)
                {
                    max_vol_constraint = dgdx;
                }  

                if (max_stress_sensitivity < min_stress_sens)
                {
                    min_stress_sens = max_stress_sensitivity;
                }

                if (max_stress_sensitivity > max_stress_sens)
                {
                    max_stress_sens = max_stress_sensitivity;
                }

            }
            
            //MIN COMPLIANCE: Fill needed vectors with the values for min compliance
            if (strcmp( optimization_type , "min_compliance" ) == 0)
            {
                for( ModelPart::ElementsContainerType::iterator element_i = mrModelPart.ElementsBegin(); element_i!= mrModelPart.ElementsEnd(); element_i++ )
                {	
                    
                    double xval = element_i->GetValue(X_PHYS_OLD);

                    //Value of the objective function sensitivity
                    double dfdx = (element_i->GetValue(DCDX));

                    // Value of the constraint function sensitivity
                    double dgdx = (element_i->GetValue(DVDX));

                    
                    double Xmin = 0;
                    double Xmax = 1;
                    vol_summ = vol_summ + xval;
                    x[iteration]= xval;
                    df[iteration]= dfdx;
                    dg[iteration] = dgdx;
                    xmax[iteration] = Xmax;
                    xmin[iteration] = Xmin;
                    iteration = iteration + 1;                        
                }
            }

            //MAX COMPLIANCE: Fill needed vectors with the values for max compliance
            if (strcmp( optimization_type , "max_compliance" ) == 0)
            {
                for( ModelPart::ElementsContainerType::iterator element_i = mrModelPart.ElementsBegin(); element_i!= mrModelPart.ElementsEnd(); element_i++ )
                {	
                    
                    double xval = element_i->GetValue(X_PHYS_OLD);

                    //Value of the objective function sensitivity
                    double dfdx = (element_i->GetValue(DCDX_COMPLIANT));

                    // Value of the constraint function sensitivity
                    double dgdx = (element_i->GetValue(DVDX));

                    double initial_size = element_i->GetValue(INITIAL_ELEMENT_SIZE);

                    if ( iteration < 2)
                        max_mean_stress = element_i->GetValue(MAX_MEAN_STRESS);

                    double max_stress_sensitivity = element_i->GetValue(YOUNGS_MODULUS_SENSITIVITY);

                    solid_void = element_i->GetValue(SOLID_VOID);
                    if (solid_void > 0 )
                        xval = 1.0;
                    
                    if ( std::isnan(dfdx))
                    {
                        KRATOS_ERROR << "The nan comes from the Element in MMA Object sens: " <<dfdx <<std::endl;
                    }

                    if ( std::isnan(dgdx))
                    {
                        KRATOS_ERROR << "The nan comes from the Element in MMA vol constraint: " <<dgdx <<std::endl;
                    }

                    if ( std::isnan(max_stress_sensitivity))
                    {
                        KRATOS_ERROR << "The nan comes from the Element in MMA stress sens: " <<max_stress_sensitivity <<std::endl;
                    }
                    
                    double Xmin = 0;
                    double Xmax = 1;
                    vol_summ = vol_summ + initial_size*xval;
                    x[iteration]= xval;
                    df[iteration]= dfdx/max_obj_func;
                    dg[iteration*mm +0] =  initial_size/(domain_size*volfrac);
                    dg[iteration*mm +1] = max_stress_sensitivity/50;
                    xmax[iteration] = std::min(Xmax, xval+move);
                    xmin[iteration] = std::max(Xmin, xval-move);
                    iteration = iteration + 1;
                    //KRATOS_INFO("[TopOpt]") << "  Sensitivität objective: "<<dfdx/max_obj_func <<"! Sens. constraint: "<< initial_size/max_vol_constraint<<"! Sens. constr. stress: "<< max_stress_sensitivity<< std::endl;         
                    
                }
            }

 			// Initialize MMA

			g[0] = 0;
			vol_frac_iteration = vol_summ;
			g[0] = ((vol_frac_iteration/(domain_size*volfrac))-1);
            g[1] = 0;
            g[1] = (max_mean_stress/50)-1;
            KRATOS_INFO("[TopOpt]") << "  constraint value: "<<g[1] << " max stress sens: "<< max_stress_sens<< " and minimum sens: " <<min_stress_sens<<std::endl;
			
			// Update the design variables using MMA

            mpMMA->Update(x,df,g,dg,xmin,xmax);


			int jiter = 0;
            double Void = 0;

			for(ModelPart::ElementsContainerType::iterator elem_i = mrModelPart.ElementsBegin();
					elem_i!=mrModelPart.ElementsEnd(); elem_i++)
				{
                Void = elem_i->GetValue(SOLID_VOID);
                if (Void == 1)
                {
                    elem_i->SetValue(X_PHYS, 1);
                    jiter= jiter +1;
                }
                else
                {
				    elem_i->SetValue(X_PHYS, x[jiter]);
				    jiter= jiter +1;
                }
				}



			// Printing of results
			clock_t end = clock();
			KRATOS_INFO("[TopOpt]") << "  Updating of values performed               [ spent time =  " << double(end - begin) / CLOCKS_PER_SEC << " ] " << std::endl;
		}
		else 
		{

			KRATOS_ERROR << "No valid optimization_algorithm selected for the simulation. Selected one: " << update_type << std::endl;
		}



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
        return "MMAAlgorithm";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "MMAAlgorithm";
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

    ModelPart& mrModelPart;
    MMASolver* mpMMA = nullptr;

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
    //MMAAlgorithm& operator=(MMAAlgorithm const& rOther);

    /// Copy constructor.
    //MMAAlgorithm(MMAAlgorithm const& rOther);


    ///@}

}; // Class MMAAlgorithm

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif	/* defined(KRATOS_MMA_ALGORITHM_H_INCLUDED */

