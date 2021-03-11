//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author Julio Marti
//

#if !defined(KRATOS_PARTICLES_UTILITIES_INCLUDED )
#define  KRATOS_PARTICLES_UTILITIES_INCLUDED

#define PRESSURE_ON_EULERIAN_MESH
#define USE_FEW_PARTICLES

// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "utilities/geometry_utilities.h"
#include "geometries/tetrahedra_3d_4.h"
#include "pfem_melting_application.h"
#include "spatial_containers/spatial_containers.h"
#include "utilities/timer.h"
#include "processes/node_erase_process.h"
#include "utilities/binbased_fast_point_locator.h"
#include "includes/deprecated_variables.h"

#include <boost/timer.hpp>
#include "utilities/timer.h"

#ifdef _OPENMP
#include "omp.h"
#endif

namespace Kratos
{
  
  template< class T, std::size_t dim >
    class DistanceCalculator1
    {
    public:
      
      double operator()(T const& p1, T const& p2)
      {
	double dist = 0.0;
        for (std::size_t i = 0; i < dim; i++)
	  {
            double tmp = p1[i] - p2[i];
            dist += tmp*tmp;
	  }
        return dist; //square distance because it is easier to work without the square root//
      }
      
    };
  
  template<std::size_t TDim> class FaceHeatFlux
    {
    public:
      KRATOS_CLASS_POINTER_DEFINITION(FaceHeatFlux<TDim>);
      
      
      void FaceHeatFluxDistribution(ModelPart & rLagrangianModelPart, double x, double y, double z, double radius, double face_heat_flux)
      {
        KRATOS_TRY
	  
	  //defintions for spatial search
	  typedef Node < 3 > PointType;
        typedef Node < 3 > ::Pointer PointTypePointer;
        typedef std::vector<PointType::Pointer> PointVector;
        typedef std::vector<PointType::Pointer>::iterator PointIterator;
        typedef std::vector<double> DistanceVector;
        typedef std::vector<double>::iterator DistanceIterator;
	
        //creating an auxiliary list for the new nodes
        PointVector list_of_nodes;

        //*************
        // Bucket types
        typedef Bucket< TDim, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > BucketType;
	
        typedef Tree< KDTreePartition<BucketType> > tree; //Kdtree;
	
	
        //starting calculating time of construction of the kdtree
        boost::timer kdtree_construction;
	
        for (ModelPart::NodesContainerType::iterator node_it = rLagrangianModelPart.NodesBegin();
	     node_it != rLagrangianModelPart.NodesEnd(); ++node_it)
	  {
            PointTypePointer pnode = *(node_it.base());
	    
            //putting the nodes of the destination_model part in an auxiliary list
            list_of_nodes.push_back(pnode);
	  }
	
        std::cout << "kdt constructin time " << kdtree_construction.elapsed() << std::endl;
	
        //create a spatial database with the list of new nodes
        unsigned int bucket_size = 20;
        tree nodes_tree(list_of_nodes.begin(), list_of_nodes.end(), bucket_size);
	
        //work arrays
        Node < 3 > work_point(0, 0.0, 0.0, 0.0);
        unsigned int MaximumNumberOfResults = 10000;
        PointVector Results(MaximumNumberOfResults);
        DistanceVector SquaredResultsDistances(MaximumNumberOfResults);
	
	
        //if (rEulerianModelPart.NodesBegin()->SolutionStepsDataHas(NODAL_H) == false)
	//  KRATOS_THROW_ERROR(std::logic_error, "Add  ----NODAL_H---- variable!!!!!! ERROR", "");
        
        double sigma = 0.0;
        if (TDim == 2)
	  sigma = 10.0 / (7.0 * 3.1415926);
        else
	  sigma = 1.0 / 3.1415926;
	
        work_point.X() = x;
	work_point.Y() = y;
	work_point.Z() = z;

	//double radius = 1.5 * node_it->FastGetSolutionStepValue(NODAL_H);
		
        //find all of the new nodes within the radius
        int number_of_points_in_radius;
		
        //look between the new nodes which of them is inside the radius of the circumscribed cyrcle
        number_of_points_in_radius = nodes_tree.SearchInRadius(work_point, radius, Results.begin(), SquaredResultsDistances.begin(), MaximumNumberOfResults);
		
        if (number_of_points_in_radius > 0)
	{
		    
	//double& temperature = (node_it)->FastGetSolutionStepValue(TEMPERATURE);
		    
		double temperature_aux = 0.0;
		    
		double tot_weight = 0.0;
		double maximunweight = SPHCubicKernel(sigma, 0.0, radius);    
		for (int k = 0; k < number_of_points_in_radius; k++)
		{
			double distance = sqrt(*(SquaredResultsDistances.begin() + k));
	
			double weight = SPHCubicKernel(sigma, distance, radius);
			
			PointIterator it_found = Results.begin() + k;
			
			if((*it_found)->FastGetSolutionStepValue(IS_FREE_SURFACE)==1) //MATERIAL_VARIABLE
			  {
			
			    double& aux= (*it_found)->FastGetSolutionStepValue(FACE_HEAT_FLUX);
                            aux =  face_heat_flux * weight / maximunweight; 
			}
		}	     
	  }
        KRATOS_CATCH("")
	  }
      
 
    private:
      
      inline double SPHCubicKernel(const double sigma, const double r, const double hmax)
      {
        double h_half = 0.5 * hmax;
        const double s = r / h_half;
        const double coeff = sigma / pow(h_half, static_cast<int>(TDim));
	
        if (s <= 1.0)
	  return coeff * (1.0 - 1.5 * s * s + 0.75 * s * s * s);
        else if (s <= 2.0)
	  return 0.25 * coeff * pow(2.0 - s, 3);
        else
            return 0.0;
      }
      


    };

} // namespace Kratos.

#endif // KRATOS_LAGRANGIAN_PARTICLES_UTILITIES_INCLUDED  defined


