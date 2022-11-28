//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes
#include "interpolative_mapper_base.h"
#include "custom_mappers/nearest_neighbor_mapper.h"
#include "custom_mappers/nearest_element_mapper.h"
#include "custom_mappers/barycentric_mapper.h"
#include "utilities/geometrical_projection_utilities.h"
#include "utilities/tessellation_utilities/delaunator_utilities.h"
#include "utilities/auxiliar_model_part_utilities.h"
#include "utilities/parallel_utilities.h"

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

/**
 * @brief Type of meta mapper considered
 */
enum class MetaMapperType
{
    NEAREST_NEIGHBOR,
    NEAREST_ELEMENT,
    BARYCENTRIC
};

///@}
///@name  Functions
///@{

namespace ModelPartUtility
{
/**
 * @brief This function returns the corresponding model part considering the parameters
 * @param rOriginModelPart The model part to be considered
 * @param ThisParameters The configuration parameters
*/
ModelPart& GetOriginModelPart(
    ModelPart& rOriginModelPart,
    Parameters ThisParameters = Parameters(R"({})")
    )
{
    // We get the model part name
    const std::string origin_2d_sub_model_part_name = ThisParameters.Has("origin_2d_sub_model_part_name") ? ThisParameters["origin_2d_sub_model_part_name"].GetString() : "";

    // We check if the submodelpart exists
    if (origin_2d_sub_model_part_name != "") {
        // We check if the submodelpart is the modelpart
        if (rOriginModelPart.Name() == origin_2d_sub_model_part_name) {
            return rOriginModelPart;
        } else { // We retrieve the orginal modelpart if not
            KRATOS_ERROR_IF_NOT(rOriginModelPart.HasSubModelPart(origin_2d_sub_model_part_name)) << "Submodelpart " << rOriginModelPart.FullName() << "." << origin_2d_sub_model_part_name << " does not exist" << std::endl;
            return rOriginModelPart.GetSubModelPart(origin_2d_sub_model_part_name);
        }
    } else { // We return the origin model part if not defined
        return rOriginModelPart;
    }
}
} // namespace Utility

///@}
///@name Kratos Classes
///@{

// Definition of projection variables
struct ProjectionVariables
{
    ProjectionVariables(array_1d<double, 3>& rNormal, Point& rPoint) : reference_normal(rNormal), reference_point(rPoint) {};
    array_1d<double, 3> reference_normal;
    Point reference_point;
    double distance;
    array_1d<double, 3> projected_point_coordinates;
};

/**
 * @ingroup MapingApplication
 * @class Projection3D2DMapper
 * @brief This mapper simplifies the mapping between two model parts thanks to the projection over a reference plane
 * @author Vicente Mataix Ferrandiz
 */
template<class TSparseSpace, class TDenseSpace, class TMapperBackend>
class Projection3D2DMapper
    : public InterpolativeMapperBase<TSparseSpace, TDenseSpace, TMapperBackend>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Projection3D2DMapper
    KRATOS_CLASS_POINTER_DEFINITION(Projection3D2DMapper);

    /// BaseType definitions
    typedef InterpolativeMapperBase<TSparseSpace, TDenseSpace, TMapperBackend> BaseType;
    typedef Kratos::unique_ptr<BaseType> BaseMapperUniquePointerType;
    typedef typename BaseType::TMappingMatrixType TMappingMatrixType;
    typedef typename BaseType::MapperUniquePointerType MapperUniquePointerType;

    /// Interface definitions
    typedef typename TMapperBackend::InterfaceCommunicatorType InterfaceCommunicatorType;
    typedef typename InterfaceCommunicator::MapperInterfaceInfoUniquePointerType MapperInterfaceInfoUniquePointerType;

    /// Other mappers definition
    typedef NearestNeighborMapper<TSparseSpace, TDenseSpace, TMapperBackend> NearestNeighborMapperType;
    typedef NearestElementMapper<TSparseSpace, TDenseSpace, TMapperBackend>   NearestElementMapperType;
    typedef BarycentricMapper<TSparseSpace, TDenseSpace, TMapperBackend>         BarycentricMapperType;

    /// Geometric definitions
    typedef Node<3> NodeType;
    typedef Geometry<NodeType> GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    // Default constructor, needed for registration
    Projection3D2DMapper(
        ModelPart& rModelPartOrigin,
        ModelPart& rModelPartDestination
        ) : BaseType(rModelPartOrigin, rModelPartDestination),
            mrInputModelPartOrigin(rModelPartOrigin),
            mrInputModelPartDestination(rModelPartDestination)
    {
    }

    // Constructor with settings
    Projection3D2DMapper(
        ModelPart& rModelPartOrigin,
        ModelPart& rModelPartDestination,
        Parameters JsonParameters
        ) : BaseType(ModelPartUtility::GetOriginModelPart(rModelPartOrigin, JsonParameters), rModelPartDestination, JsonParameters),
            mrInputModelPartOrigin(rModelPartOrigin),
            mrInputModelPartDestination(rModelPartDestination)
    {
        KRATOS_TRY;

        // Validate input
        this->ValidateInput();

        // Copying parameters
        mCopiedParameters = JsonParameters.Clone();

        // Type of metamapper considered
        GetBaseMapperType();
        const bool is_geometric_based_mapper = (mMetaMapperType == MetaMapperType::NEAREST_ELEMENT || mMetaMapperType == MetaMapperType::NEAREST_ELEMENT) ? true :  false;

        // 2D model parts name if any
        const std::string origin_2d_sub_model_part_name = mCopiedParameters["origin_2d_sub_model_part_name"].GetString();
        mIs2DOrigin = (origin_2d_sub_model_part_name != "") ? true : false;

        KRATOS_ERROR_IF(!mIs2DOrigin && is_geometric_based_mapper && mrInputModelPartOrigin.IsDistributed()) << "ERROR:: Geometric based mappers (\"nearest_element\"  or \"barycentric\") only work in MPI when passing a 2D model part as origin modelpart, as DelaunatorUtilities work only in serial" << std::endl;

        // We retrieve the values of interest
        if (!mIs2DOrigin) {
            noalias(mNormalPlane) = mCopiedParameters["normal_plane"].GetVector();
            noalias(mPointPlane.Coordinates()) = mCopiedParameters["reference_plane_coordinates"].GetVector();
        } else {
            const auto& r_origin_sub_model_part = this->GetOriginModelPart();
            GeometryType::Pointer p_geometry = nullptr;
            int entity = 0;
            if (r_origin_sub_model_part.NumberOfElements() > 0) {
                const auto first_element = r_origin_sub_model_part.ElementsBegin();
                p_geometry = first_element->pGetGeometry();
            } else if (r_origin_sub_model_part.NumberOfConditions() > 0) {
                const auto first_condition = r_origin_sub_model_part.ConditionsBegin();
                p_geometry = first_condition->pGetGeometry();
                entity = 1;
            }

            // MPI data
            int aux_partition_entity = -1;
            const auto& r_communicator = r_origin_sub_model_part.GetCommunicator();
            const auto& r_data_communicator = r_communicator.GetDataCommunicator();
            const int mpi_rank = r_data_communicator.Rank();
            const int mpi_size = r_data_communicator.Size();

            // Get maximum rank
            if (p_geometry != nullptr) aux_partition_entity = mpi_rank;

            // Define the send tag
            const int tag_send_index = 1;
            const int tag_send_normal = 2;
            const int tag_send_point = 3;

            // Determine root's rank
            const int root_rank = 0;

            // Getting the partition with entities index
            int partition_entity = r_data_communicator.Max(aux_partition_entity, root_rank);

            // Transfer to other partitions
            if (mrInputModelPartOrigin.IsDistributed()) {
                if (mpi_rank == root_rank) {
                    for (int i_rank = 1; i_rank < mpi_size; ++i_rank) {
                        r_data_communicator.Send(partition_entity, i_rank, tag_send_index);
                    }
                } else {
                    r_data_communicator.Recv(partition_entity, root_rank, tag_send_index);
                }
            }
 
            // Getting from parameters if not elements or conditions
            if (partition_entity != mpi_rank) {
                if (!mrInputModelPartOrigin.IsDistributed()) {
                    KRATOS_ERROR_IF_NOT(mMetaMapperType == MetaMapperType::NEAREST_NEIGHBOR) << "The mapper \"nearest_element\"  or \"barycentric\" cannot be used without elements or conditions" << std::endl; // NOTE: In the MPI case if the number of elements and conditions is zero also the number of nodes is zero and therefore not projection is needed, so we can proceed even with default values of normals and reference points
                    noalias(mNormalPlane) = mCopiedParameters["normal_plane"].GetVector();
                    noalias(mPointPlane.Coordinates()) = mCopiedParameters["reference_plane_coordinates"].GetVector();
                } else { // Now transfer the normal plane and the point plane between partitions
                    // The partitions that receive
                    r_data_communicator.Recv(mNormalPlane, partition_entity, tag_send_normal);
                    r_data_communicator.Recv(mPointPlane.Coordinates(), partition_entity, tag_send_point);
                }
            } else {
                GeometryType::CoordinatesArrayType aux_coords;
                noalias(mPointPlane.Coordinates()) = p_geometry->Center();
                p_geometry->PointLocalCoordinates(aux_coords, mPointPlane);
                noalias(mNormalPlane) = p_geometry->UnitNormal(aux_coords);

                // Doing a check that all normals are aligned
                std::size_t check_normal;
                const double numerical_limit = std::numeric_limits<double>::epsilon() * 1.0e4;
                struct normal_check {
                    normal_check(array_1d<double, 3>& rNormal) : reference_normal(rNormal) {};
                    array_1d<double, 3> reference_normal;
                    GeometryType::CoordinatesArrayType aux_coords;
                };
                if (entity == 0) {
                    check_normal = block_for_each<SumReduction<std::size_t>>(r_origin_sub_model_part.Elements(), normal_check(mNormalPlane), [&numerical_limit](auto& r_elem, normal_check& nc) {
                        auto& r_geom = r_elem.GetGeometry();
                        r_geom.PointLocalCoordinates(nc.aux_coords, r_geom.Center());
                        const auto normal = r_geom.UnitNormal(nc.aux_coords);
                        if (norm_2(normal - nc.reference_normal) > numerical_limit) { return 1; } else { return 0; }
                    });
                } else {
                    check_normal = block_for_each<SumReduction<std::size_t>>(r_origin_sub_model_part.Conditions(), normal_check(mNormalPlane), [&numerical_limit](auto& r_cond, normal_check& nc) {
                        auto& r_geom = r_cond.GetGeometry();
                        r_geom.PointLocalCoordinates(nc.aux_coords, r_geom.Center());
                        const auto normal = r_geom.UnitNormal(nc.aux_coords);
                        if (norm_2(normal - nc.reference_normal) > numerical_limit) { return 1; } else { return 0; }
                    });
                }
                KRATOS_ERROR_IF_NOT(check_normal == 0) << "The 2D reference model part has not consistent normals. Please check that is properly aligned" << std::endl;

                // The partition that sends
                if (mrInputModelPartOrigin.IsDistributed()) {
                    const auto& r_point_coordinates = mPointPlane.Coordinates();
                    for (int i_rank = 0; i_rank < mpi_size; ++i_rank) {
                        if (i_rank != partition_entity) {
                            r_data_communicator.Send(mNormalPlane, i_rank, tag_send_normal);
                            r_data_communicator.Send(r_point_coordinates, i_rank, tag_send_point);
                        }
                    }
                }
            }
        }

        // Cleaning the parameters
        mCopiedParameters.RemoveValue("normal_plane");
        mCopiedParameters.RemoveValue("reference_plane_coordinates");
        mCopiedParameters.RemoveValue("base_mapper");
        mCopiedParameters.RemoveValue("origin_2d_sub_model_part_name");
        mCopiedParameters.RemoveValue("destination_2d_sub_model_part_name");

        // Move mesh
        MoveModelParts();

        // Initializing the base mapper
        CreateBaseMapper();

        // Unmove mesh
        UnMoveModelParts();

        // Calling initialize
        this->Initialize();

        // Now we copy the mapping matrix
        BaseType::mpMappingMatrix = Kratos::make_unique<TMappingMatrixType>(mpBaseMapper->GetMappingMatrix());

        // Remove the created model parts for inverted mapping conflict
        auto& r_origin_model = mrInputModelPartOrigin.GetModel();
        if (r_origin_model.HasModelPart("projected_origin_modelpart")) {
            r_origin_model.DeleteModelPart("projected_origin_modelpart");
        }

        KRATOS_CATCH("");
    }

    /// Destructor.
    ~Projection3D2DMapper() override = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Updates the mapping-system after the geometry/mesh has changed
     * After changes in the topology (e.g. remeshing or sliding interfaces)
     * the relations for the mapping have to be recomputed. This means that
     * the search has to be conducted again and the mapping-system has to be
     * rebuilt, hence this is expensive
     * @param MappingOptions flags used to specify how the update has to be done
     * @param SearchRadius search radius used for the search
     */
    void UpdateInterface(
        Kratos::Flags MappingOptions,
        double SearchRadius
        ) override
    {
        KRATOS_TRY;
        
        // Move mesh
        MoveModelParts();

        // Initializing the base mapper
        CreateBaseMapper();

        // Update interface base mapper
        mpBaseMapper->UpdateInterface(MappingOptions, SearchRadius);

        // Unmove mesh
        UnMoveModelParts();

        // Calling initialize
        BaseType::UpdateInterface(MappingOptions, SearchRadius);

        // Now we copy the mapping matrix
        BaseType::mpMappingMatrix = Kratos::make_unique<TMappingMatrixType>(mpBaseMapper->GetMappingMatrix());

        // Remove the created model parts for inverted mapping conflict
        auto& r_origin_model = mrInputModelPartOrigin.GetModel();
        if (r_origin_model.HasModelPart("projected_origin_modelpart")) {
            r_origin_model.DeleteModelPart("projected_origin_modelpart");
        }

        KRATOS_CATCH("");
    }

    /**
    * @brief Cloning the Mapper
    * returns a clone of the current Mapper
    * pure virtual, has to be implemented in every derived mapper,
    * used in the creation of the Mappers
    * @see MapperFactory
    */
    MapperUniquePointerType Clone(
        ModelPart& rModelPartOrigin,
        ModelPart& rModelPartDestination,
        Parameters JsonParameters
        ) const override
    {
        KRATOS_TRY;

        return Kratos::make_unique<Projection3D2DMapper<TSparseSpace, TDenseSpace, TMapperBackend>>(rModelPartOrigin, rModelPartDestination, JsonParameters);

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    int AreMeshesConforming() const override
    {
        KRATOS_WARNING_ONCE("Mapper") << "Developer-warning: \"AreMeshesConforming\" is deprecated and will be removed in the future" << std::endl;
        return BaseType::mMeshesAreConforming;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "Projection3D2DMapper";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Projection3D2DMapper";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        BaseType::PrintData(rOStream);
    }

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

    /**
     * @brief This function creates the inverted mapping parameters if they are required to be differemt from the forward mapping parameters
     * @details This function has to be implemented in the derived classes in case the inverted mapping parameters are required to be different from the forward mapping parameters
     * @return The inverted mapping parameters
     */
    Parameters GetInvertedMappingParameters(Parameters ForwardMappingParameters) override
    {
        KRATOS_TRY;

        // Copy the parameters
        Parameters inverted_parameters = ForwardMappingParameters.Clone();

        // Invserse 2D submodelparts mapping parameters
        const std::string origin_2d_sub_model_part_name = inverted_parameters["origin_2d_sub_model_part_name"].GetString();
        const std::string destination_2d_sub_model_part_name = inverted_parameters["destination_2d_sub_model_part_name"].GetString();
        inverted_parameters["origin_2d_sub_model_part_name"].SetString(destination_2d_sub_model_part_name);
        inverted_parameters["destination_2d_sub_model_part_name"].SetString(origin_2d_sub_model_part_name);
        
        return inverted_parameters;

        KRATOS_CATCH("");
    }

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

    ModelPart& mrInputModelPartOrigin;          /// The unaltered origin model part (as it comes from the input)
    ModelPart& mrInputModelPartDestination;     /// The unaltered destination model part (as it comes from the input)
    BaseMapperUniquePointerType mpBaseMapper;   /// Pointer to the base mapper
    array_1d<double, 3> mNormalPlane;           /// The normal defining the plane to project
    Point mPointPlane;                          /// The coordinates of the plane to project
    Parameters mCopiedParameters;               /// The copied parameters. We copy the parameters to avoid conflicts in inverse mapping
    MetaMapperType mMetaMapperType;             /// The meta mapper type
    bool mIs2DOrigin;                           /// If the origin model part is 2D

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Getting the base mapper type
     */
    void GetBaseMapperType()
    {
        KRATOS_TRY;

        // Detect the mapper type
        const std::string mapper_name = mCopiedParameters["base_mapper"].GetString();
        if (mapper_name == "nearest_neighbor") {
            mMetaMapperType = MetaMapperType::NEAREST_NEIGHBOR;
        } else if (mapper_name == "nearest_element") {
            mMetaMapperType = MetaMapperType::NEAREST_ELEMENT;
        } else if (mapper_name == "barycentric") {
            mMetaMapperType = MetaMapperType::BARYCENTRIC;
        } else {
            KRATOS_ERROR << "Not defined mapper type " << mapper_name << std::endl;
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Create the base mapper
     */
    void CreateBaseMapper()
    {
        KRATOS_TRY;

        // The origin model         
        auto& r_origin_model = mrInputModelPartOrigin.GetModel();

        // Getting model parts
        auto& r_projected_origin_modelpart = r_origin_model.HasModelPart("projected_origin_modelpart") ? r_origin_model.GetModelPart("projected_origin_modelpart") : this->GetOriginModelPart();
        auto& r_projected_destination_modelpart = mrInputModelPartDestination;

        // Creating the base mapper
        if (mMetaMapperType == MetaMapperType::NEAREST_NEIGHBOR) {
            if (mCopiedParameters.Has("interpolation_type")) mCopiedParameters.RemoveValue("interpolation_type");
            if (mCopiedParameters.Has("local_coord_tolerance")) mCopiedParameters.RemoveValue("local_coord_tolerance");
            mpBaseMapper = Kratos::make_unique<NearestNeighborMapperType>(r_projected_origin_modelpart, r_projected_destination_modelpart, mCopiedParameters);
        } else if (mMetaMapperType == MetaMapperType::NEAREST_ELEMENT) {
            if (mCopiedParameters.Has("interpolation_type")) mCopiedParameters.RemoveValue("interpolation_type");
            mpBaseMapper = Kratos::make_unique<NearestElementMapperType>(r_projected_origin_modelpart, r_projected_destination_modelpart, mCopiedParameters);
        } else if (mMetaMapperType == MetaMapperType::NEAREST_ELEMENT) {
            mpBaseMapper = Kratos::make_unique<BarycentricMapperType>(r_projected_origin_modelpart, r_projected_destination_modelpart, mCopiedParameters);
        } else {
            KRATOS_ERROR << "ERROR:: Mapper " << mCopiedParameters["base_mapper"].GetString() << " is not available as base mapper for projection" << std::endl;
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Moves the model parts
     */
    void MoveModelParts()
    {
        KRATOS_TRY;

        // Projected origin model part. If the origin model part is not 2D we project it
        if (!mIs2DOrigin) {
            const bool is_geometric_based_mapper = (mMetaMapperType == MetaMapperType::NEAREST_ELEMENT || mMetaMapperType == MetaMapperType::NEAREST_ELEMENT) ? true :  false;
            if(is_geometric_based_mapper) {
                // Origin model
                auto& r_origin_model = mrInputModelPartOrigin.GetModel();

                auto& r_projected_origin_modelpart = r_origin_model.CreateModelPart("projected_origin_modelpart");

                // Iterate over the existing nodes
                double distance;
                array_1d<double, 3> projected_point_coordinates;
                for (auto& r_node : mrInputModelPartOrigin.Nodes()) {
                    noalias(projected_point_coordinates) = GeometricalProjectionUtilities::FastProject(mPointPlane, r_node, mNormalPlane, distance).Coordinates();
                    r_projected_origin_modelpart.CreateNewNode(r_node.Id(), projected_point_coordinates[0], projected_point_coordinates[1], projected_point_coordinates[2]);
                }

                // We generate "geometries" to be able to interpolate
                DelaunatorUtilities::CreateTriangleMeshFromNodes(r_projected_origin_modelpart);
            } else { // Move mesh
                MapperUtilities::SaveCurrentConfiguration(mrInputModelPartOrigin);

                // Iterate over the existing nodes
                block_for_each(mrInputModelPartOrigin.Nodes(), ProjectionVariables(mNormalPlane, mPointPlane), [&](auto& r_node, ProjectionVariables& p) {
                    noalias(p.projected_point_coordinates) = GeometricalProjectionUtilities::FastProject(p.reference_point, r_node, p.reference_normal, p.distance).Coordinates();
                    noalias(r_node.Coordinates()) = p.projected_point_coordinates;
                });
            }
        } 

        // Projected destination model part
        MapperUtilities::SaveCurrentConfiguration(mrInputModelPartDestination);

        // Iterate over the existing nodes
        block_for_each(mrInputModelPartDestination.Nodes(), ProjectionVariables(mNormalPlane, mPointPlane), [&](auto& r_node, ProjectionVariables& p) {
            noalias(p.projected_point_coordinates) = GeometricalProjectionUtilities::FastProject(p.reference_point, r_node, p.reference_normal, p.distance).Coordinates();
            noalias(r_node.Coordinates()) = p.projected_point_coordinates;
        });

        KRATOS_CATCH("");
    }

    /**
     * @brief Unmoves the model parts
     */
    void UnMoveModelParts()
    {
        KRATOS_TRY; 

        const bool is_geometric_based_mapper = (mMetaMapperType == MetaMapperType::NEAREST_ELEMENT || mMetaMapperType == MetaMapperType::NEAREST_ELEMENT) ? true :  false;
        if (!mIs2DOrigin && !is_geometric_based_mapper) {
            MapperUtilities::RestoreCurrentConfiguration(mrInputModelPartOrigin);
        }
        MapperUtilities::RestoreCurrentConfiguration(mrInputModelPartDestination);

        KRATOS_CATCH("");
    }

    /**
     * @brief This function origin model part (for inverse mapping)
     * @return The origin model part
     */
    ModelPart& GetOriginModelPartForInverseMapping() override
    {
        return mrInputModelPartDestination;
    }

    /**
     * @brief This function destination model part (for inverse mapping)
     * @return The destination model part
     */
    ModelPart& GetDestinationModelPartForInverseMapping() override
    {
        return mrInputModelPartOrigin;
    }

    void CreateMapperLocalSystems(
        const Communicator& rModelPartCommunicator,
        std::vector<Kratos::unique_ptr<MapperLocalSystem>>& rLocalSystems
        ) override
    {
        // Calling base mapper method. But not sure if this must be changed
        AccessorInterpolativeMapperBase<TMapperBackend>::CreateMapperLocalSystems(*mpBaseMapper, rModelPartCommunicator, rLocalSystems);
    }

    MapperInterfaceInfoUniquePointerType GetMapperInterfaceInfo() const override
    {
        return AccessorInterpolativeMapperBase<TMapperBackend>::GetMapperInterfaceInfo(*mpBaseMapper);
    }

    Parameters GetMapperDefaultSettings() const override
    {
        return Parameters( R"({
            "search_settings"                    : {},
            "echo_level"                         : 0,
            "interpolation_type"                 : "unspecified",
            "origin_2d_sub_model_part_name"      : "",
            "destination_2d_sub_model_part_name" : "",
            "local_coord_tolerance"              : 0.25,
            "use_initial_configuration"          : false,
            "print_pairing_status_to_file"       : false,
            "pairing_status_file_path"           : "",
            "base_mapper"                        : "nearest_neighbor",
            "normal_plane"                       : [0.0,0.0,1.0],
            "reference_plane_coordinates"        : [0.0,0.0,0.0]
        })");
    }

    ///@}

}; // Class Projection3D2DMapper

///@} addtogroup block
}  // namespace Kratos.