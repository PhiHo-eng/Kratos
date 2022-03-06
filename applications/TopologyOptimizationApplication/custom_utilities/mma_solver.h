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

    void UpdateDensitiesUsingMMAAlgorithm( char update_type[], double volfrac, double greyscale , double OptItr , double qmax)
    {
        KRATOS_TRY;

        if ( strcmp( update_type , "mma_algorithm" ) == 0 )
		{
			clock_t begin = clock();
			KRATOS_INFO("[TopOpt]") << "  Method of Moving Asymthodes (MMA) chosen to solve   the optimization problem" << std::endl;


			// Create object of updating function
			int nn = mrModelPart.NumberOfElements();
			int mm = 1; /// constraints

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

			//get the information from the element and save them in the new vectors
			for( ModelPart::ElementsContainerType::iterator element_i = mrModelPart.ElementsBegin(); element_i!= mrModelPart.ElementsEnd(); element_i++ )
			{	
				
				double xval = element_i->GetValue(X_PHYS);

				//Value of the objective function sensitivity
				double dfdx = (element_i->GetValue(DCDX_COMPLIANT));

				// Value of the constraint function sensitivity
				double dgdx = (element_i->GetValue(DVDX));


				double youngs_modulus = element_i->GetProperties()[YOUNGS_MODULUS_0];

				
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


 			// Initialize MMA

			g[0] = 0;
			vol_frac_iteration = vol_summ;
			g[0] = (vol_frac_iteration -  volfrac*nn);

			double Xminn = 0;
			double Xmaxx= 1;
			double movlim = 0.5;
	

			for (int iEl = 0; iEl < nn; iEl++) 
			{
				xmax[iEl] = std::min(Xmaxx, x[iEl] + movlim);
				xmin[iEl] = std::max(Xminn, x[iEl] - movlim); 
			}
			
			// Update the design variables using MMA

            mpMMA->Update(x,df,g,dg,xmin,xmax);


			int jiter = 0;

			for(ModelPart::ElementsContainerType::iterator elem_i = mrModelPart.ElementsBegin();
					elem_i!=mrModelPart.ElementsEnd(); elem_i++)
				{
				elem_i->SetValue(X_PHYS, x[jiter]);
				jiter= jiter +1;
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

