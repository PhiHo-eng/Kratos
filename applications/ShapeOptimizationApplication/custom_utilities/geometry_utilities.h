// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef GEOMETRY_UTILITIES_H
#define GEOMETRY_UTILITIES_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <functional>

#include <pybind11/pybind11.h>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/key_hash.h"
#include "shape_optimization_application.h"
#include "utilities/variable_utils.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "includes/global_pointer_variables.h"
#include "utilities/element_size_calculator.h"
#include "utilities/atomic_utilities.h"

#include "spatial_containers/spatial_containers.h"

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

class GeometryUtilities
{
public:
    ///@name Type Definitions
    ///@{

    // For better reading
    typedef array_1d<double,3> array_3d;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;
    typedef ModelPart::ElementType::GeometryType GeometryType;
    typedef std::size_t SizeType;

    /// Pointer definition of GeometryUtilities
    KRATOS_CLASS_POINTER_DEFINITION(GeometryUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    GeometryUtilities( ModelPart& modelPart )
        : mrModelPart( modelPart )
    {
    }

    /// Destructor.
    virtual ~GeometryUtilities()
    {
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // --------------------------------------------------------------------------
    void ComputeUnitSurfaceNormals()
    {
        KRATOS_TRY;

        const unsigned int domain_size = mrModelPart.GetProcessInfo().GetValue(DOMAIN_SIZE);
        KRATOS_ERROR_IF(mrModelPart.NumberOfConditions() == 0) <<
            "> Normal calculation requires surface or line conditions to be defined!" << std::endl;
        KRATOS_ERROR_IF((domain_size == 3 && mrModelPart.ConditionsBegin()->GetGeometry().size() == 2)) <<
            "> Normal calculation of 2-noded conditions in 3D domains is not possible!" << std::endl;
        CalculateAreaNormalsFromConditions();
        CalculateUnitNormals();

        KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------
    void ProjectNodalVariableOnUnitSurfaceNormals( const Variable<array_3d> &rNodalVariable )
    {
        KRATOS_TRY;

        ProjectNodalVariableOnDirection(rNodalVariable, NORMALIZED_SURFACE_NORMAL);

        KRATOS_CATCH("");
    }

    void ProjectNodalVariableOnDirection( const Variable<array_3d> &rNodalVariable, const Variable<array_3d> &rDirectionVariable)
    {
        KRATOS_TRY;

        for (ModelPart::NodeIterator node_i = mrModelPart.NodesBegin(); node_i != mrModelPart.NodesEnd(); ++node_i)
        {
            array_3d &nodal_variable = node_i->FastGetSolutionStepValue(rNodalVariable);
            array_3d &node_normal = node_i->FastGetSolutionStepValue(rDirectionVariable);

            const double magnitude = inner_prod(nodal_variable, node_normal);
            noalias(nodal_variable) = magnitude * node_normal;
        }

        KRATOS_CATCH("");
    }

    void ProjectNodalVariableOnTangentPlane( const Variable<array_3d> &rNodalVariable, const Variable<array_3d> &rPlaneNormalVariable)
    {
        KRATOS_TRY;

        for (ModelPart::NodeIterator node_i = mrModelPart.NodesBegin(); node_i != mrModelPart.NodesEnd(); ++node_i)
        {
            array_3d &nodal_variable = node_i->FastGetSolutionStepValue(rNodalVariable);
            array_3d &node_normal = node_i->FastGetSolutionStepValue(rPlaneNormalVariable);

            const double magnitude = inner_prod(nodal_variable, node_normal);
            nodal_variable -= magnitude * node_normal;
        }

        KRATOS_CATCH("");
    }
    // --------------------------------------------------------------------------
    void ExtractBoundaryNodes( std::string const& rBoundarySubModelPartName )
    {
    	KRATOS_TRY;

        ModelPart& r_boundary_model_part = mrModelPart.GetSubModelPart(rBoundarySubModelPartName);

        KRATOS_ERROR_IF(r_boundary_model_part.Nodes().size() != 0) << "ExtractBoundaryNodes: The boundary model part already has nodes!" << std::endl;

    	// Some type-definitions
        typedef std::unordered_map<vector<unsigned int>, unsigned int, KeyHasherRange<vector<unsigned int>>, KeyComparorRange<vector<unsigned int>> > hashmap;

    	// Create map to ask for number of faces for the given set of node ids representing one face in the model part
    	hashmap n_boundaries_map;

        unsigned int domain_size = static_cast<unsigned int>(mrModelPart.GetProcessInfo()[DOMAIN_SIZE]);

    	// Fill map that counts number of faces for given set of nodes
    	for (auto& elem_i : mrModelPart.Elements())
    	{
            KRATOS_ERROR_IF(elem_i.GetGeometry().Dimension() < domain_size) << "ExtractBoundaryNodes: This function does only work"
                <<" for solid elements in 3D and surface elements in 2D!" << std::endl;

            Element::GeometryType::GeometriesArrayType boundaries = elem_i.GetGeometry().GenerateBoundariesEntities();

            for(unsigned int boundary=0; boundary<boundaries.size(); boundary++)
            {
                // Create vector that stores all node is of current face
                DenseVector<unsigned int> ids(boundaries[boundary].size());

                // Store node ids
                for(unsigned int i=0; i<boundaries[boundary].size(); i++)
                    ids[i] = boundaries[boundary][i].Id();

                //*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
                std::sort(ids.begin(), ids.end());

                // Fill the map
                n_boundaries_map[ids] += 1;
            }
    	}

    	// Vector to store all nodes on surface. Node ids may be listed several times
    	std::vector<std::size_t> temp_boundary_node_ids;

    	// Add surface nodes to sub-model part
    	for(auto it=n_boundaries_map.begin(); it!=n_boundaries_map.end(); it++)
    	{
    		// If given node set represents face that is not overlapping with a face of another element, add it as skin element
    		if(it->second == 1)
    		{
    			for(unsigned int i=0; i<it->first.size(); i++)
    				temp_boundary_node_ids.push_back(it->first[i]);
    		}
    	}

    	// Add nodes and remove double entries
    	r_boundary_model_part.AddNodes(temp_boundary_node_ids);

    	KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------
    void ExtractEdgeNodes( std::string const& rEdgeSubModelPartName ) {
        KRATOS_TRY;

        KRATOS_ERROR_IF(mrModelPart.Elements().size() == 0) << "ExtractEdgeNodes: No elements defined. Automatic edge detection will not find any edge nodes!" << std::endl;

        ModelPart& r_edge_model_part = mrModelPart.GetSubModelPart(rEdgeSubModelPartName);

        KRATOS_ERROR_IF(r_edge_model_part.Nodes().size() != 0) << "ExtractEdgeNodes: The edge model part already has nodes!" << std::endl;

        for (auto& r_node_i : mrModelPart.Nodes()) {
            auto& r_node_i_neighbours = r_node_i.GetValue(NEIGHBOUR_NODES);  // does not work with const
            for(const auto& r_node_j : r_node_i_neighbours) {
                auto& r_element_neighbours = r_node_i.GetValue(NEIGHBOUR_ELEMENTS);  // does not work with const
                int count = 0;
                for(const auto& r_elem_k : r_element_neighbours) {
                    const auto& r_element_geometry = r_elem_k.GetGeometry();
                    for(const auto& r_node_l : r_element_geometry) {
                        if (r_node_l.Id() == r_node_j.Id()) count ++;
                    }
                }
                if (count < 2){
                    r_edge_model_part.AddNode(&r_node_i);
                    break;
                }
            }
        }
        KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------
    void ComputeDistancesToBoundingModelPart(
        ModelPart& rBoundingModelPart,
        pybind11::list& rSignedDistances,
        pybind11::list& rDirections )
    {
        KRATOS_TRY;

        typedef Node < 3 > NodeType;
        typedef NodeType::Pointer NodeTypePointer;
        typedef std::vector<NodeType::Pointer> NodeVector;
        typedef std::vector<NodeType::Pointer>::iterator NodeIterator;
        typedef std::vector<double>::iterator DoubleVectorIterator;
        typedef Bucket< 3, NodeType, NodeVector, NodeTypePointer, NodeIterator, DoubleVectorIterator > BucketType;
        typedef Tree< KDTreePartition<BucketType> > KDTree;

        KRATOS_ERROR_IF(rBoundingModelPart.NumberOfElements() != 0) <<
            "ComputeDistancesToBoundingModelPart: Model part must only contain conditions!" << std::endl;

        NodeVector all_bounding_nodes;
        all_bounding_nodes.reserve(rBoundingModelPart.Nodes().size());
        for (ModelPart::NodesContainerType::iterator node_it = rBoundingModelPart.NodesBegin(); node_it != rBoundingModelPart.NodesEnd(); ++node_it)
        {
            all_bounding_nodes.push_back(*(node_it.base()));
        }
        const size_t bucket_size = 100;
        KDTree search_tree(all_bounding_nodes.begin(), all_bounding_nodes.end(), bucket_size);

        GeometryUtilities(rBoundingModelPart).ComputeUnitSurfaceNormals();

        for (auto& r_node : mrModelPart.Nodes()){

            double distance;
            NodeTypePointer p_neighbor = search_tree.SearchNearestPoint(r_node, distance);

            const array_3d delta = r_node.Coordinates() - p_neighbor->Coordinates();
            const array_3d& bounding_normal = p_neighbor->FastGetSolutionStepValue(NORMALIZED_SURFACE_NORMAL);
            const double projected_length = inner_prod(delta, bounding_normal);

            rSignedDistances.append(projected_length);

            rDirections.append(bounding_normal[0]);
            rDirections.append(bounding_normal[1]);
            rDirections.append(bounding_normal[2]);
        }

        KRATOS_CATCH("");
    }


    template <class TContainerType>
    double CalculateLength(TContainerType& rContainer)
    {
        double length = block_for_each<SumReduction<double>>(rContainer, [&](typename TContainerType::value_type& rEntity){
            return rEntity.GetGeometry().Length();
        });

        return length;
    }

// --------------------------------------------------------------------------
    double ComputeVolume()
    {
        KRATOS_TRY

        const double volume = block_for_each<SumReduction<double>>(
            mrModelPart.Elements(), [&](const ModelPart::ElementType& rElement) {
                return rElement.GetGeometry().DomainSize();
            });

        return mrModelPart.GetCommunicator().GetDataCommunicator().SumAll(volume);

        KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------
    void ComputeVolumeShapeDerivatives(
        const Variable<array_3d>& rDerivativeVariable)
    {
        KRATOS_TRY

        using CUInt = const unsigned int;

        KRATOS_ERROR_IF_NOT(mrModelPart.HasNodalSolutionStepVariable(rDerivativeVariable))
            << rDerivativeVariable.Name() << " not found in solution step variables list in "
            << mrModelPart.FullName() << ".\n";

        VariableUtils().SetHistoricalVariableToZero(rDerivativeVariable, mrModelPart.Nodes());

        const int domain_size = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];
        const SizeType number_of_nodes = GetNumberOfElementalNodesInModelPart(mrModelPart);

        std::function<double(CUInt, CUInt, const GeometryType&)> volume_derivative_method;

        if (domain_size == 2) {
            if (number_of_nodes == 3) {
                volume_derivative_method = [](CUInt NodeIndex, CUInt DirectionIndex,  const GeometryType& rGeometry) {
                    return 2.0 * ElementSizeCalculator<2, 3>::AverageElementSize(rGeometry) * ElementSizeCalculator<2, 3>::AverageElementSizeDerivative(NodeIndex, DirectionIndex, rGeometry);
                };
            } else if (number_of_nodes == 4) {
                volume_derivative_method = [](CUInt NodeIndex, CUInt DirectionIndex,  const GeometryType& rGeometry) {
                    return 2.0 * ElementSizeCalculator<2, 4>::AverageElementSize(rGeometry) * ElementSizeCalculator<2, 4>::AverageElementSizeDerivative(NodeIndex, DirectionIndex, rGeometry);
                };
            } else {
                KRATOS_ERROR
                    << "Unsupported geometry type. Only geometries with 3 and "
                       "4 nodes are supported in 2D. [ number_of_nodes = "
                    << number_of_nodes << " ].\n";
            }
        } else if (domain_size == 3) {
            if (number_of_nodes == 4) {
                volume_derivative_method = [](CUInt NodeIndex, CUInt DirectionIndex,  const GeometryType& rGeometry) {
                    return 3.0 * std::pow(ElementSizeCalculator<3, 4>::AverageElementSize(rGeometry), 2) * ElementSizeCalculator<3, 4>::AverageElementSizeDerivative(NodeIndex, DirectionIndex, rGeometry);
                };
            } else if (number_of_nodes == 6) {
                volume_derivative_method = [](CUInt NodeIndex, CUInt DirectionIndex,  const GeometryType& rGeometry) {
                    return 3.0 * std::pow(ElementSizeCalculator<3, 6>::AverageElementSize(rGeometry), 2) * ElementSizeCalculator<3, 6>::AverageElementSizeDerivative(NodeIndex, DirectionIndex, rGeometry);
                };
            } else if (number_of_nodes == 8) {
                volume_derivative_method = [](CUInt NodeIndex, CUInt DirectionIndex,  const GeometryType& rGeometry) {
                    return 3.0 * std::pow(ElementSizeCalculator<3, 8>::AverageElementSize(rGeometry), 2) * ElementSizeCalculator<3, 8>::AverageElementSizeDerivative(NodeIndex, DirectionIndex, rGeometry);
                };
            } else {
                KRATOS_ERROR
                    << "Unsupported geometry type. Only geometries with 4, 6 and "
                       "8 nodes are supported in 3D. [ number_of_nodes = "
                    << number_of_nodes << " ].\n";
            }
        } else {
            KRATOS_ERROR << "Unsupported domain size. Only 2D and 3D is "
                            "supported [ domain_size = "
                         << domain_size << " ].\n";
        }

        block_for_each(mrModelPart.Elements(), [&](ModelPart::ElementType& rElement){
            auto& r_geometry = rElement.GetGeometry();
            for (SizeType c = 0; c < r_geometry.PointsNumber(); ++c) {
                auto& r_derivative_value = r_geometry[c].FastGetSolutionStepValue(rDerivativeVariable);

                for (int k = 0; k < domain_size; ++k) {
                    const double derivative_value = volume_derivative_method(c, k, r_geometry);
                    AtomicAdd(r_derivative_value[k], derivative_value);
                }
            }
        });

        mrModelPart.GetCommunicator().AssembleCurrentData(rDerivativeVariable);

        KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------

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
        return "GeometryUtilities";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "GeometryUtilities";
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

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    // --------------------------------------------------------------------------
    void CalculateAreaNormalsFromConditions()
    {
        KRATOS_TRY

        //resetting the normals
        VariableUtils().SetHistoricalVariableToZero(NORMAL, mrModelPart.Nodes());

        //calculating the normals and summing up at nodes
        const array_1d<double,3> local_coords = ZeroVector(3);
        block_for_each(mrModelPart.Conditions(), [&](Condition& rCond)
        {
            auto& r_geometry = rCond.GetGeometry();

            const array_1d<double,3> normal = r_geometry.Normal(local_coords);
            const double coeff = 1.0/r_geometry.size();

            for(auto& node_i : r_geometry)
            {
                node_i.SetLock();
                noalias(node_i.FastGetSolutionStepValue(NORMAL)) += coeff * normal;
                node_i.UnSetLock();
            }
        });

        KRATOS_CATCH("")
    }

    // --------------------------------------------------------------------------
    void CalculateUnitNormals()
    {
        for (auto& node_i : mrModelPart.Nodes())
        {
            const array_1d<double,3>& area_normal = node_i.FastGetSolutionStepValue(NORMAL);
            array_3d& normalized_normal = node_i.FastGetSolutionStepValue(NORMALIZED_SURFACE_NORMAL);

            const double norm2 = norm_2(area_normal);
            KRATOS_ERROR_IF(norm2<1e-10) << "CalculateUnitNormals: Norm2 of normal for node "
                << node_i.Id() << " is < 1e-10!" << std::endl;

            noalias(normalized_normal) = area_normal/norm2;
        }
    }

    // --------------------------------------------------------------------------
    SizeType GetNumberOfElementalNodesInModelPart(const ModelPart& rModelPart) const
    {
        KRATOS_TRY

        const int local_geometry_type = block_for_each<MaxReduction<int>>(mrModelPart.Elements(), [&](const ModelPart::ElementType& rElement) {
            return static_cast<int>(rElement.GetGeometry().GetGeometryType());
        });
        const int global_geometry_type = mrModelPart.GetCommunicator().GetDataCommunicator().MaxAll(local_geometry_type);

        block_for_each(mrModelPart.Elements(), [&](const ModelPart::ElementType& rElement) {
            KRATOS_ERROR_IF(static_cast<int>(rElement.GetGeometry().GetGeometryType()) != global_geometry_type)
                << "Mismatching geometries found in "
                << mrModelPart.Name() << ". All geometries should be of same type. [ First element's geometry's type = "
                << global_geometry_type << ", current_element_id = " << rElement.Id()
                << ", current element's geometry's type = "
                << static_cast<int>(rElement.GetGeometry().GetGeometryType()) << " ].\n";
        });

        const SizeType local_geometry_number_of_nodes = block_for_each<MaxReduction<SizeType>>(mrModelPart.Elements(), [&](const ModelPart::ElementType& rElement) {
            return rElement.GetGeometry().PointsNumber();
        });
        return mrModelPart.GetCommunicator().GetDataCommunicator().MaxAll(local_geometry_number_of_nodes);

        KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------


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
//      GeometryUtilities& operator=(GeometryUtilities const& rOther);

    /// Copy constructor.
//      GeometryUtilities(GeometryUtilities const& rOther);


    ///@}

}; // Class GeometryUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // GEOMETRY_UTILITIES_H
