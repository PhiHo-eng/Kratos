//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#pragma once

// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes

// Project includes
#include "containers/edge_based_data_structure.h"
#include "includes/define.h"
#include "includes/gid_io.h"
#include "processes/process.h"
#include "utilities/atomic_utilities.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

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
/**takes a model part full of SIMPLICIAL ELEMENTS (triangles and tetras) and convects a level set distance
* on the top of it
*/
template<unsigned int TDim>
class FluxCorrectedTransportConvectionProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    using SizeType = std::size_t;

    using IndexType = std::size_t;

    ///@}
    ///@name Pointer Definitions
    ///@{

    /// Pointer definition of FluxCorrectedTransportConvectionProcess
    KRATOS_CLASS_POINTER_DEFINITION(FluxCorrectedTransportConvectionProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    FluxCorrectedTransportConvectionProcess(
        Model& rModel,
        Parameters ThisParameters)
        : Process()
    {
        // Validate the common settings as well as the element formulation specific ones
        ThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());

        // Reference to the model part to which the convection is to be performed
        mpModelPart = &(rModel.GetModelPart(ThisParameters["model_part_name"].GetString()));

        // Set the edge-based data structure pointer
        mpEdgeDataStructure = Kratos::make_unique<EdgeBasedDataStructure<TDim>>();

        // Assign all the required member variables
        mEchoLevel = ThisParameters["echo_level"].GetInt();
        mMaxAllowedCFL = ThisParameters["max_CFL"].GetDouble();
        mMaxAllowedDt = ThisParameters["max_delta_time"].GetDouble();
        mpConvectedVar = &KratosComponents<Variable<double>>::Get(ThisParameters["convected_variable_name"].GetString());
        mpConvectionVar = &KratosComponents<Variable<array_1d<double, 3>>>::Get(ThisParameters["convection_variable_name"].GetString());
    }

    /// Copy constructor.
    FluxCorrectedTransportConvectionProcess(FluxCorrectedTransportConvectionProcess const& rOther) = delete;

    /// Destructor.
    ~FluxCorrectedTransportConvectionProcess() = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    FluxCorrectedTransportConvectionProcess &operator=(FluxCorrectedTransportConvectionProcess const &rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    void ExecuteInitialize() override
    {
        // Fill the edge-based data structure
        mpEdgeDataStructure->CalculateEdgeDataStructure(*mpModelPart);

        // Auxiliary flag to avoid resetting the edge data structure at each Execute call
        mPerformInitialize = false;

        // Allocate auxiliary arrays
        // Note that we give the size such that the index matches the id of the corresponding node
        // Also note that this implies having a "fake" 0-id entry in the first position of the arrays
        mAuxSize = mpModelPart->NumberOfNodes() + 1;
        mResidual.resize(mAuxSize);
        mSolution.resize(mAuxSize);
        mSolutionOld.resize(mAuxSize);
        mLowOrderUpdate.resize(mAuxSize);
        mHighOrderUpdate.resize(mAuxSize);
        mConvectionValues.resize(mAuxSize);
    }

    void Execute() override
    {
        KRATOS_TRY;

        // Check if the edge data structure is already set
        if (mPerformInitialize) {
            this->ExecuteInitialize();
        }

        // Evaluate current time (n) and target time (n+1)
        // Note that it is assumed that the CloneTimeStep of the model part has been already done
        const auto& r_process_info = mpModelPart->GetProcessInfo();
        const auto &r_prev_process_info = r_process_info.GetPreviousTimeStepInfo();
        const double time = r_process_info[TIME];
        double prev_time = r_prev_process_info[TIME];

        // Initialize data vectors
        const SizeType n_nodes = mpModelPart->NumberOfNodes();
        IndexPartition<IndexType>(n_nodes).for_each([this](IndexType iNode){
            const auto it_node = mpModelPart->NodesBegin() + iNode;
            const IndexType i_node_id = it_node->Id();
            mSolutionOld[i_node_id] = it_node->GetSolutionStepValue(*mpConvectedVar,1); //solution at n
            mConvectionValues[i_node_id] = it_node->GetSolutionStepValue(*mpConvectionVar,1); //convective velocity at n
        });

        GidIO<> gid_io_convection("/home/rzorrilla/Desktop/FluxCorrectedTransportConvectionProcess2D", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
        gid_io_convection.InitializeMesh(0);
        gid_io_convection.WriteMesh(mpModelPart->GetMesh());
        gid_io_convection.FinalizeMesh();
        gid_io_convection.InitializeResults(0, mpModelPart->GetMesh());
        gid_io_convection.WriteNodalResults(DISTANCE, mpModelPart->Nodes(), 0, 0);
        gid_io_convection.WriteNodalResults(VELOCITY, mpModelPart->Nodes(), 0, 0);

        // Substepping loop according to max CFL
        IndexType step = 1;
        while (time - prev_time > 1.0e-12) {
            // Evaluate current substep maximum allowable delta time
            const double dt = this->CalculateSubStepDeltaTime(prev_time, time);

            // Solve current substep
            KRATOS_INFO_IF("FluxCorrectedTransportConvectionProcess", mEchoLevel > 0) << "Substep " << step << " - \u0394t " << dt << std::endl;
            this->SolveSubStep(dt);

            // TODO:Remove after debugging
            IndexPartition<IndexType>(mpModelPart->NumberOfNodes()).for_each([this](IndexType iNode){
                auto it_node = mpModelPart->NodesBegin() + iNode;
                it_node->FastGetSolutionStepValue(*mpConvectedVar) = mSolution[it_node->Id()];
            });
            gid_io_convection.WriteNodalResults(DISTANCE, mpModelPart->Nodes(), step, 0);
            gid_io_convection.WriteNodalResults(VELOCITY, mpModelPart->Nodes(), step, 0);

            // Advance in time
            step++;
            prev_time += dt;
        }

        gid_io_convection.FinalizeResults();

        // Set final solution in the model part database
        IndexPartition<IndexType>(mpModelPart->NumberOfNodes()).for_each([this](IndexType iNode){
            auto it_node = mpModelPart->NodesBegin() + iNode;
            it_node->FastGetSolutionStepValue(*mpConvectedVar) = mSolution[it_node->Id()];
        });

        KRATOS_CATCH("")
    }

    void Clear() override
    {
        mResidual.clear();
        mSolution.clear();
        mSolutionOld.clear();
        mLowOrderUpdate.clear();
        mHighOrderUpdate.clear();
        mConvectionValues.clear();
        mpEdgeDataStructure->Clear();
    }

    const Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"({
            "echo_level" : 0,
            "model_part_name" : "",
            "convected_variable_name" : "DISTANCE",
            "convection_variable_name" : "VELOCITY",
            "max_CFL" : 1.0,
            "max_delta_time" : 1.0
        })");

        return default_parameters;
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
    std::string Info() const override
    {
        return "FluxCorrectedTransportConvectionProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "FluxCorrectedTransportConvectionProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ModelPart* mpModelPart; /// Pointer to the background mesh model part

    SizeType mAuxSize; /// Size of the internal vector containers that equals the number of nodes plus one

    SizeType mEchoLevel; /// Level of information that is output

    double mMaxAllowedDt; /// Maximum allowed time step for the substepping

    double mMaxAllowedCFL; /// Maximum allowed CFL number for the substepping

    bool mPerformInitialize = true; /// Flag to indicate if the ExecuteInitialize is required to be done

    const Variable<double>* mpConvectedVar = nullptr; /// Pointer to the convected variable (e.g. DISTANCE)

    const Variable<array_1d<double,3>>* mpConvectionVar = nullptr; /// Pointer to the convection variable (e.g. VELOCITY)

    std::vector<double> mResidual; /// Auxiliary vector to assemble the edge residuals

    std::vector<double> mSolution; /// Auxiliary vector containing the current solution

    std::vector<double> mSolutionOld; /// Auxiliary vector containing the previous step solution

    std::vector<double> mLowOrderUpdate; /// Auxiliary vector containing the low order solution update

    std::vector<double> mHighOrderUpdate; /// Auxiliary vector containing the high order solution update

    std::vector<array_1d<double,3>> mConvectionValues; /// Auxiliary vector to store the convection variable values (e.g. VELOCITY)

    typename EdgeBasedDataStructure<TDim>::UniquePointer mpEdgeDataStructure; /// Pointer to the edge-based data structure

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    //TODO: Do it with the edge projected velocity
    double CalculateSubStepDeltaTime(
        const double PreviousTime,
        const double TargetTime)
    {
        // Get edge data structure containers
        const auto& r_row_indices = mpEdgeDataStructure->GetRowIndices();
        const auto& r_col_indices = mpEdgeDataStructure->GetColIndices();

        // Calculate the CFL number at each edge to get the minimum allowed time step
        // Note that we don't loop the last entry of the row container as it is NNZ
        SizeType n_nodes = r_row_indices.size() - 1;
        const double min_dt = IndexPartition<IndexType>(n_nodes).for_each<MinReduction<double>>([&](IndexType iNode){
            double dt_ij = std::numeric_limits<double>::max();
            const auto it_row = r_row_indices.begin() + iNode;
            const IndexType i_col_index = *it_row;
            const SizeType n_cols = *(it_row+1) - i_col_index;
            // Check that there are CSR columns (i.e. that current node involves an edge)
            if (n_cols != 0) {
                // i-node nodal data
                const double i_vel_norm = norm_2(mConvectionValues[iNode]);

                // j-node nodal loop (i.e. loop ij-edges)
                const auto i_col_begin = r_col_indices.begin() + i_col_index;
                for (IndexType j_node = 0; j_node < n_cols; ++j_node) {
                    // j-node nodal data
                    IndexType j_node_id = *(i_col_begin + j_node);
                    const double j_vel_norm = norm_2(mConvectionValues[j_node_id]);

                    // Get ij-edge length from CSR data structure
                    const auto& r_ij_edge_data = mpEdgeDataStructure->GetEdgeData(iNode, j_node_id);
                    const double ij_length = r_ij_edge_data.GetLength();

                    // Calculate the max time step from the velocities at both edge ends and get the minimum
                    dt_ij = CalculateEdgeLocalDeltaTime(ij_length, i_vel_norm, j_vel_norm);
                }
            }

            return dt_ij;
        });

        // Check that current time step doesn't exceed the target time
        // Use maximum allowable time step if it doesn't exceed the target time
        // Otherwise the difference between current time and target one is used (potentially in last substep)
        double dt = PreviousTime + min_dt > TargetTime ? TargetTime - PreviousTime : min_dt;

        // Synchronize among MPI nodes
        mpModelPart->GetCommunicator().GetDataCommunicator().MinAll(dt);

        // Return maximum allowable time step
        return dt;
    }

    double CalculateEdgeLocalDeltaTime(
        const double Length,
        const double NormVelI,
        const double NormVelJ)
    {
        const double dx_CFL = Length * mMaxAllowedCFL;
        const double dt_i = NormVelI > 1.0e-12 ? dx_CFL / NormVelI : mMaxAllowedDt;
        const double dt_j = NormVelJ > 1.0e-12 ? dx_CFL / NormVelJ : mMaxAllowedDt;
        return std::min(dt_i, dt_j);
    }

    void SolveSubStep(const double DeltaTime)
    {
        // Initialize containers for the current substep solve
        IndexPartition<IndexType>(mAuxSize).for_each([this](IndexType i){
            mResidual[i] = 0.0;
            mLowOrderUpdate[i] = 0.0;
            mHighOrderUpdate[i] = 0.0;
        });

        // Calculate residual
        CalculateResidual(DeltaTime);

        // Calculate the low order solution
        CalculateLowOrderUpdate(DeltaTime);

        // Calculate the high order solution update
        CalculateHighOrderSolutionUpdate(DeltaTime);

        // Calculate the antidiffusive edge contributions
        CalculateAntidiffusiveEdgeContributions(DeltaTime);

        // Add the low order update
        // TODO: Check that this needs to be done for the limiter or not
        IndexPartition<IndexType>(mAuxSize).for_each([this](IndexType i){
            mSolution[i] = mSolutionOld[i] + mLowOrderUpdate[i];
        });

        // Evaluate limiter
        EvaluateLimiter(DeltaTime);

        // Do the current step solution update
        IndexPartition<IndexType>(mAuxSize).for_each([this](IndexType i){
            mSolutionOld[i] = mSolution[i];
        });
    }

    void CalculateResidual(const double DeltaTime)
    {
        // Get edge data structure containers
        const auto &r_row_indices = mpEdgeDataStructure->GetRowIndices();
        const auto &r_col_indices = mpEdgeDataStructure->GetColIndices();
        SizeType aux_n_rows = r_row_indices.size() - 1; // Note that the last entry of the row container is the NNZ

        // Aux TLS container
        struct AuxTLS
        {
            array_1d<double, TDim> F_i;
            array_1d<double, TDim> F_j;
            array_1d<double, TDim> d_ij;
            array_1d<double, TDim> b_ij;
            array_1d<double, TDim> F_ij_num;
            array_1d<double, 3> vel_ij_half;
        };

        // Off-diagonal (edge) contributions assembly
        IndexPartition<IndexType>(aux_n_rows).for_each(AuxTLS(), [&](IndexType iRow, AuxTLS& rTLS){
            // Get current row (node) storage data
            const auto it_row = r_row_indices.begin() + iRow;
            const IndexType i_col_index = *it_row;
            const SizeType n_cols = *(it_row+1) - i_col_index;

            // Check that there are CSR columns (i.e. that current node involves an edge)
            if (n_cols != 0) {
                // Get TLS data
                auto& F_i = rTLS.F_i;
                auto& F_j = rTLS.F_j;
                auto& d_ij = rTLS.d_ij;
                auto& b_ij = rTLS.b_ij;
                auto& F_ij_num = rTLS.F_ij_num;
                auto& vel_ij_half = rTLS.vel_ij_half;

                // i-node nodal data
                const double u_i = mSolutionOld[iRow];
                const auto& r_i_vel = mConvectionValues[iRow];

                // i-node convective flux calculation
                for (IndexType d = 0; d < TDim; ++d) {
                    F_i[d] = u_i * r_i_vel[d];
                }

                // j-node nodal loop (i.e. loop ij-edges)
                const auto i_col_begin = r_col_indices.begin() + i_col_index;
                for (IndexType j_node = 0; j_node < n_cols; ++j_node) {
                    // j-node nodal data
                    IndexType j_node_id = *(i_col_begin + j_node);
                    const double u_j = mSolutionOld[j_node_id];
                    const auto& r_j_vel = mConvectionValues[j_node_id];

                    // j-node convective flux calculation
                    for (IndexType d = 0; d < TDim; ++d) {
                        F_j[d] = u_j * r_j_vel[d];
                    }

                    // Get ij-edge operators data from CSR data structure
                    const auto& r_ij_edge_data = mpEdgeDataStructure->GetEdgeData(iRow, j_node_id);
                    const auto& r_Ni_DNj = r_ij_edge_data.GetOffDiagonalConvective();
                    const auto& r_DNi_Nj = r_ij_edge_data.GetOffDiagonalConvectiveTranspose();
                    d_ij = 0.5 * (r_Ni_DNj - r_DNi_Nj);
                    double D_ij = 0.0;
                    for (IndexType d = 0; d < TDim; ++d) {
                        D_ij += std::pow(d_ij[d],2);
                    }

                    // Calculate numerical flux "upwind" contribution at de edge midpoint
                    // Note that this numerical flux corresponds to the Lax-Wendroff scheme
                    const double u_ij_half = 0.5 * (u_i + u_j) - 0.5 * DeltaTime * inner_prod(-d_ij, F_i - F_j) / D_ij;
                    noalias(vel_ij_half) = 0.5 * (r_i_vel + r_j_vel);
                    for (IndexType d = 0; d < TDim; ++d) {
                        F_ij_num[d] = 2.0 * (vel_ij_half[d] * u_ij_half);
                    }

                    // Calculate convection volume residual contributions
                    double res_edge_i = inner_prod(-d_ij, F_ij_num);
                    double res_edge_j = inner_prod(d_ij, F_ij_num);

                    // If current edge belogs to a boundary, add the corresponding convection boundary integrals
                    if (r_ij_edge_data.IsBoundary()) {
                        // Get ij-edge boundary operators from CSR data structure
                        const auto &r_N_N_normal = r_ij_edge_data.GetConvectiveBoundary();
                        b_ij = 0.5 * r_N_N_normal;

                        // Add boundary contribution to the residual
                        const double res_edge_bd = inner_prod(b_ij, F_ij_num);
                        res_edge_i -= res_edge_bd;
                        res_edge_j += res_edge_bd;
                    }

                    // Atomic additions
                    AtomicAdd(mResidual[iRow], res_edge_i);
                    AtomicAdd(mResidual[j_node_id], res_edge_j);
                }
            }
        });

        // Add the diagonal boundary term coming from the convective flux split
        IndexPartition<IndexType>(mpModelPart->NumberOfNodes()).for_each(array_1d<double,TDim>(), [&](IndexType iNode, array_1d<double,TDim>& rTLS){
            // Get nodal data
            const auto it_node = mpModelPart->NodesBegin() + iNode;
            const IndexType i_node_id = it_node->Id();

            // Add the diagonal contribution to the residual
            rTLS = mpEdgeDataStructure->GetBoundaryMassMatrixDiagonal(i_node_id);
            double aux_res = 0.0;
            const double u_i = mSolutionOld[i_node_id];
            const auto &r_i_vel = mConvectionValues[i_node_id];
            for (IndexType d = 0; d < TDim; ++d) {
                aux_res += u_i * r_i_vel[d] * rTLS[d];
            }
            mResidual[i_node_id] -= aux_res;
        });
    }

    void CalculateLowOrderUpdate(const double DeltaTime)
    {
        // Get edge data structure containers
        const auto &r_row_indices = mpEdgeDataStructure->GetRowIndices();
        const auto &r_col_indices = mpEdgeDataStructure->GetColIndices();
        SizeType aux_n_rows = r_row_indices.size() - 1; // Note that the last entry of the row container is the NNZ

        // Initialize low order update container to the Taylor-Galerkin residual
        IndexPartition<IndexType>(mAuxSize).for_each([this](IndexType i){
            mLowOrderUpdate[i] = mResidual[i];
        });

        // Add the diffusion to the Taylor-Galerkin already computed residual
        IndexPartition<IndexType>(aux_n_rows).for_each([&](IndexType iRow){
            // Get current row (node) storage data
            const auto it_row = r_row_indices.begin() + iRow;
            const IndexType i_col_index = *it_row;
            const SizeType n_cols = *(it_row+1) - i_col_index;

            // Check that there are CSR columns (i.e. that current node involves an edge)
            if (n_cols != 0) {
                // i-node nodal data
                const double u_i = mSolutionOld[iRow];
                const auto& r_i_vel = mConvectionValues[iRow];

                // j-node nodal loop (i.e. loop ij-edges)
                const auto i_col_begin = r_col_indices.begin() + i_col_index;
                for (IndexType j_node = 0; j_node < n_cols; ++j_node) {
                    // j-node nodal data
                    IndexType j_node_id = *(i_col_begin + j_node);
                    const double u_j = mSolutionOld[j_node_id];
                    const auto& r_j_vel = mConvectionValues[j_node_id];

                    // Add low order scheme diffusion
                    const auto& r_ij_edge_data = mpEdgeDataStructure->GetEdgeData(iRow, j_node_id);
                    const double local_dt = CalculateEdgeLocalDeltaTime(r_ij_edge_data.GetLength(), norm_2(r_i_vel), norm_2(r_j_vel));
                    const double c_tau = 1.0 / local_dt;
                    const double Mc_ij = r_ij_edge_data.GetOffDiagonalConsistentMass();
                    const double ij_low_order_diff =  c_tau * Mc_ij;
                    const double res_edge_i = ij_low_order_diff * (u_j - u_i);
                    const double res_edge_j = ij_low_order_diff * (u_i - u_j);

                    // Atomic additions
                    AtomicAdd(mLowOrderUpdate[iRow], res_edge_i);
                    AtomicAdd(mLowOrderUpdate[j_node_id], res_edge_j);
                }
            }
        });

        // Do the explicit lumped mass matrix solve to obtain the low order update
        IndexPartition<IndexType>(mpModelPart->NumberOfNodes()).for_each([&](IndexType iNode){
            // Get nodal data
            const auto it_node = mpModelPart->NodesBegin() + iNode;
            const IndexType i_node_id = it_node->Id();

            // Do the explicit lumped mass matrix solve
            const double M_l = mpEdgeDataStructure->GetLumpedMassMatrixDiagonal(i_node_id);
            const double solve_coeff = DeltaTime / M_l;
            mLowOrderUpdate[i_node_id] *= solve_coeff;
        });
    }

    void CalculateHighOrderSolutionUpdate(const double DeltaTime)
    {
        // Get edge data structure containers
        const auto &r_row_indices = mpEdgeDataStructure->GetRowIndices();
        const auto &r_col_indices = mpEdgeDataStructure->GetColIndices();
        SizeType aux_n_rows = r_row_indices.size() - 1; // Note that the last entry of the row container is the NNZ

        // Allocate auxiliar high order residual vector
        std::vector<double> residual_high_order(mAuxSize);

        // Do the consistent mass matrix solve (iteratively with lumped mass matrix)
        const SizeType max_it = 5;
        for (IndexType i = 0; i < max_it; ++i) {

            // Initialize current iteration residual to the Galerkin one
            IndexPartition<IndexType>(mAuxSize).for_each([&](IndexType i){
                residual_high_order[i] = mResidual[i];
            });

            // Add the lumping correction
            IndexPartition<IndexType>(aux_n_rows).for_each([&](IndexType iRow){
                // Get current row (node) storage data
                const auto it_row = r_row_indices.begin() + iRow;
                const IndexType i_col_index = *it_row;
                const SizeType n_cols = *(it_row+1) - i_col_index;

                // Check that there are CSR columns (i.e. that current node involves an edge)
                if (n_cols != 0) {
                    // i-node nodal data
                    const double delta_u_h_i = mHighOrderUpdate[iRow];

                    // j-node nodal loop (i.e. loop ij-edges)
                    const auto i_col_begin = r_col_indices.begin() + i_col_index;
                    for (IndexType j_node = 0; j_node < n_cols; ++j_node) {
                        // j-node nodal data
                        IndexType j_node_id = *(i_col_begin + j_node);
                        const double delta_u_h_j = mHighOrderUpdate[j_node_id];

                        // Add previous iteration correction
                        const auto& r_ij_edge_data = mpEdgeDataStructure->GetEdgeData(iRow, j_node_id);
                        const double Mc_ij = r_ij_edge_data.GetOffDiagonalConsistentMass();
                        const double res_edge_i = Mc_ij * (delta_u_h_i - delta_u_h_j) / DeltaTime;
                        const double res_edge_j = Mc_ij * (delta_u_h_j - delta_u_h_i) / DeltaTime;

                        // Atomic additions
                        AtomicAdd(residual_high_order[iRow], res_edge_i);
                        AtomicAdd(residual_high_order[j_node_id], res_edge_j);
                    }
                }
            });

            // Solve current iteration residual
            IndexPartition<IndexType>(mpModelPart->NumberOfNodes()).for_each([&](IndexType iNode){
                // Get nodal data
                const auto it_node = mpModelPart->NodesBegin() + iNode;
                const IndexType i_node_id = it_node->Id();

                // Do the explicit lumped mass matrix solve
                const double M_l = mpEdgeDataStructure->GetLumpedMassMatrixDiagonal(i_node_id);
                const double solve_coeff = DeltaTime / M_l;
                mHighOrderUpdate[i_node_id] = solve_coeff * residual_high_order[i_node_id];
            });
        }
    }

    void CalculateAntidiffusiveEdgeContributions(const double DeltaTime)
    {
        // Get edge data structure containers
        const auto &r_row_indices = mpEdgeDataStructure->GetRowIndices();
        const auto &r_col_indices = mpEdgeDataStructure->GetColIndices();
        SizeType aux_n_rows = r_row_indices.size() - 1; // Note that the last entry of the row container is the NNZ

        // Add the diffusion to the Taylor-Galerkin already computed residual
        IndexPartition<IndexType>(aux_n_rows).for_each([&](IndexType iRow){
            // Get current row (node) storage data
            const auto it_row = r_row_indices.begin() + iRow;
            const IndexType i_col_index = *it_row;
            const SizeType n_cols = *(it_row+1) - i_col_index;

            // Check that there are CSR columns (i.e. that current node involves an edge)
            if (n_cols != 0) {
                // i-node nodal data
                const double u_i = mSolutionOld[iRow];
                const double u_h_i = mHighOrderUpdate[iRow];
                const auto& r_i_vel = mConvectionValues[iRow];

                // j-node nodal loop (i.e. loop ij-edges)
                const auto i_col_begin = r_col_indices.begin() + i_col_index;
                for (IndexType j_node = 0; j_node < n_cols; ++j_node) {
                    // j-node nodal data
                    IndexType j_node_id = *(i_col_begin + j_node);
                    const double u_j = mSolutionOld[j_node_id];
                    const double u_h_j = mHighOrderUpdate[j_node_id];
                    const auto& r_j_vel = mConvectionValues[j_node_id];

                    // Calculate and store the antidiffusive flux edge contribution
                    auto& r_ij_edge_data = mpEdgeDataStructure->GetEdgeData(iRow, j_node_id);
                    const double local_dt = CalculateEdgeLocalDeltaTime(r_ij_edge_data.GetLength(), norm_2(r_i_vel), norm_2(r_j_vel));
                    const double c_tau = DeltaTime / local_dt;
                    const double Mc_ij = r_ij_edge_data.GetOffDiagonalConsistentMass();
                    const double AEC_ij = Mc_ij * (c_tau * (u_i - u_j) + (u_h_i - u_h_j)) / DeltaTime;
                    r_ij_edge_data.SetAntidiffusiveEdgeContribution(AEC_ij);
                }
            }
        });
    }

    void EvaluateLimiter(const double DeltaTime)
    {
        // Allocate and initialize auxiliary limiter arrays
        std::vector<double> P_min(mAuxSize, 0.0);
        std::vector<double> P_max(mAuxSize, 0.0);
        std::vector<double> Q_min(mAuxSize, 0.0);
        std::vector<double> Q_max(mAuxSize, 0.0);
        std::vector<double> u_min(mAuxSize, std::numeric_limits<double>::max());
        std::vector<double> u_max(mAuxSize, std::numeric_limits<double>::lowest());

        // Get edge data structure containers
        const auto &r_row_indices = mpEdgeDataStructure->GetRowIndices();
        const auto &r_col_indices = mpEdgeDataStructure->GetColIndices();
        SizeType aux_n_rows = r_row_indices.size() - 1; // Note that the last entry of the row container is the NNZ

        // Calculate the summation of antidiffusive fluxes and bound fluxes
        IndexPartition<IndexType>(aux_n_rows).for_each([&](IndexType iRow){
            // Get current row (node) storage data
            const auto it_row = r_row_indices.begin() + iRow;
            const IndexType i_col_index = *it_row;
            const SizeType n_cols = *(it_row+1) - i_col_index;

            // Check that there are CSR columns (i.e. that current node involves an edge)
            if (n_cols != 0) {
                // i-node nodal data
                const double u_i = mSolutionOld[iRow];
                const double M_l_i = mpEdgeDataStructure->GetLumpedMassMatrixDiagonal(iRow);

                // j-node nodal loop (i.e. loop ij-edges)
                const auto i_col_begin = r_col_indices.begin() + i_col_index;
                for (IndexType j_node = 0; j_node < n_cols; ++j_node) {
                    // j-node nodal data
                    IndexType j_node_id = *(i_col_begin + j_node);
                    const double u_j = mSolutionOld[j_node_id];

                    // Get the antidiffusive flux edge contribution
                    const auto& r_ij_edge_data = mpEdgeDataStructure->GetEdgeData(iRow, j_node_id);
                    const double AEC_ij = r_ij_edge_data.GetAntidiffusiveEdgeContribution();

                    // Assemble the positive/negative antidiffusive fluxes
                    AtomicAdd(P_min[iRow], std::min(0.0, AEC_ij));
                    AtomicAdd(P_max[iRow], std::max(0.0, AEC_ij));
                    AtomicAdd(P_min[j_node_id], std::min(0.0, -AEC_ij));
                    AtomicAdd(P_max[j_node_id], std::max(0.0, -AEC_ij));

                    // Assemble the positive/negative antidiffusive bounds
                    const double aux_val = M_l_i * (u_j - u_i) / DeltaTime;
                    Q_min[iRow] = std::min(Q_min[iRow], aux_val);
                    Q_max[iRow] = std::max(Q_max[iRow], aux_val);
                    Q_min[j_node_id] = std::min(Q_min[j_node_id], -aux_val);
                    Q_max[j_node_id] = std::max(Q_max[j_node_id], -aux_val);
                }
            }
        });

        // Calculate the nodal correction factors
        std::vector<double> R_min(mAuxSize);
        std::vector<double> R_max(mAuxSize);
        IndexPartition<IndexType>(mAuxSize).for_each([&](IndexType i){
            const double P_min_i = P_min[i];
            const double P_max_i = P_max[i];
            if (P_max_i > 1.0e-12 && P_min_i < -1.0e-12) {
                R_min[i] = std::min(1.0, Q_min[i] / P_min_i);
                R_max[i] = std::min(1.0, Q_max[i] / P_max_i);
            } else {
                R_min[i] = 0.0;
                R_max[i] = 0.0;
            }
        });

        // Assembly of the final antidiffusive contribution
        std::vector<double> antidiff_assembly(mAuxSize, 0.0);
        IndexPartition<IndexType>(aux_n_rows).for_each([&](IndexType iRow){
            // Get current row (node) storage data
            const auto it_row = r_row_indices.begin() + iRow;
            const IndexType i_col_index = *it_row;
            const SizeType n_cols = *(it_row+1) - i_col_index;

            // Check that there are CSR columns (i.e. that current node involves an edge)
            if (n_cols != 0) {
                // j-node nodal loop (i.e. loop ij-edges)
                const auto i_col_begin = r_col_indices.begin() + i_col_index;
                for (IndexType j_node = 0; j_node < n_cols; ++j_node) {
                    // j-node nodal data
                    IndexType j_node_id = *(i_col_begin + j_node);

                    // Get the antidiffusive flux edge contribution
                    const auto& r_ij_edge_data = mpEdgeDataStructure->GetEdgeData(iRow, j_node_id);
                    const double AEC_ij = r_ij_edge_data.GetAntidiffusiveEdgeContribution();

                    // Set as correction factor the minimum in the edge
                    const double alpha_ij = AEC_ij > 0.0 ? std::min(R_max[iRow], R_min[j_node_id]) : std::min(R_min[iRow], R_max[j_node_id]);

                    // Assemble current edge antidiffusive contribution
                    const double aux = alpha_ij * AEC_ij;
                    AtomicAdd(antidiff_assembly[iRow], aux);
                    AtomicAdd(antidiff_assembly[j_node_id], -aux);
                }
            }
        });

        // Solve and add final antidiffusive contributions
        // Note that in here it is assumed that the low order solution is already added
        IndexPartition<IndexType>(mpModelPart->NumberOfNodes()).for_each([&](IndexType iNode){
            // Get nodal data
            const auto it_node = mpModelPart->NodesBegin() + iNode;
            const IndexType i_node_id = it_node->Id();

            // Solve antidiffusive contribution
            const double M_l = mpEdgeDataStructure->GetLumpedMassMatrixDiagonal(i_node_id);
            const double solve_coeff = DeltaTime / M_l;
            mSolution[i_node_id] += solve_coeff * antidiff_assembly[i_node_id];
        });
    }

    // void CalculateLocalBounds(
    //     std::vector<double>& rMinLocalBoundVector,
    //     std::vector<double>& rMaxLocalBoundVector)
    // {
    //     // Calculate maximum and minimum allowed values
    //     block_for_each(mpModelPart->Nodes(), [&](Node &rNode) {
    //         // Initialize values for current node (i)
    //         double u_min = std::numeric_limits<double>::max();
    //         double u_max = std::numeric_limits<double>::lowest();

    //         // i-node data
    //         const IndexType i_id = rNode.Id();
    //         const double u_i_low = mSolution[i_id];
    //         const double u_i_old = mSolutionOld[i_id];
    //         const double u_i_min = std::min(u_i_low, u_i_old);
    //         const double u_i_max = std::max(u_i_low, u_i_old);

    //         // j-node nodal loop (i.e. loop ij-edges)
    //         auto& r_neighs_i = rNode.GetValue(NEIGHBOUR_NODES);
    //         for (auto& rp_neigh : r_neighs_i) {
    //             // j-node data
    //             const IndexType j_id = rp_neigh.Id();
    //             const double u_j_low = mSolution[j_id];
    //             const double u_j_old = mSolutionOld[j_id];
    //             const double u_j_min = std::min(u_j_low, u_j_old);
    //             const double u_j_max = std::max(u_j_low, u_j_old);

    //             // Check among current ij-edge nodes
    //             const double u_ij_min = std::min(u_j_min, u_i_min);
    //             const double u_ij_max = std::max(u_j_max, u_i_max);

    //             // Check among edges sharing i-node
    //             u_min = std::min(u_min, u_ij_min);
    //             u_max = std::max(u_max, u_ij_max);
    //         }

    //         // Save maximum and minimum of all edges surrounding i-node
    //         rMinLocalBoundVector[i_id] = u_min;
    //         rMaxLocalBoundVector[i_id] = u_max;
    //     });
    // }

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}
}; // Class FluxCorrectedTransportConvectionProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// Input stream function
template<unsigned int TDim>
inline std::istream& operator >> (
    std::istream& rIStream,
    FluxCorrectedTransportConvectionProcess<TDim>& rThis);

/// Output stream function
template< unsigned int TDim>
inline std::ostream& operator << (
    std::ostream& rOStream,
    const FluxCorrectedTransportConvectionProcess<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.
