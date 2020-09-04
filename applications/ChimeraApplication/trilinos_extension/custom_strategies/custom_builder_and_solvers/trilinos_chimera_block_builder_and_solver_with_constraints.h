//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Aditya Ghantasala
//
//

#if !defined(TRILINOS_CHIMERA_BLOCK_BUILDER_AND_SOLVER_WITH_CONSTRAINTS)
#define TRILINOS_CHIMERA_BLOCK_BUILDER_AND_SOLVER_WITH_CONSTRAINTS

/* System includes */
#include <set>

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "custom_strategies/builder_and_solvers/trilinos_block_builder_and_solver.h"
#include "utilities/timer.h"

/* Trilinos includes */
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_Import.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"

#define START_TIMER(label, rank) \
    if (BaseType::mrComm.MyPID() == rank)  \
        Timer::Start(label);
#define STOP_TIMER(label, rank) \
    if (BaseType::mrComm.MyPID() == rank) \
        Timer::Stop(label);

namespace Kratos {

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

/**
 * @class TrilinosChimeraBlockBuilderAndSolver
 * @ingroup TrilinosApplication
 * @brief Current class provides an extension to the trilinos b&s with constraints
 * @details
 * @author Aditya Ghantasala
 */
template <class TSparseSpace,
          class TDenseSpace,  //= DenseSpace<double>,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
class TrilinosChimeraBlockBuilderAndSolver
    : public TrilinosBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> {
public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION(TrilinosChimeraBlockBuilderAndSolver);

    /// Definition of the base class
    typedef TrilinosBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    /// The size_t types
    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    /// Definition of the classes from the base class
    typedef typename BaseType::TSchemeType TSchemeType;
    typedef typename BaseType::TDataType TDataType;
    typedef typename BaseType::DofsArrayType DofsArrayType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;
    typedef typename BaseType::NodesArrayType NodesArrayType;
    typedef typename BaseType::ElementsArrayType ElementsArrayType;
    typedef typename BaseType::ConditionsArrayType ConditionsArrayType;
    typedef typename BaseType::ElementsContainerType ElementsContainerType;

    /// Epetra definitions
    typedef Epetra_MpiComm EpetraCommunicatorType;

    /// DoF types definition
    typedef Node<3> NodeType;
    typedef typename NodeType::DofType DofType;
    typedef DofType::Pointer DofPointerType;

    typedef Element::EquationIdVectorType EquationIdVectorType;
    typedef Element::DofsVectorType DofsVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    TrilinosChimeraBlockBuilderAndSolver(EpetraCommunicatorType& rComm,
                                  int GuessRowSize,
                                  typename TLinearSolver::Pointer pNewLinearSystemSolver)
        : BaseType(rComm, GuessRowSize, pNewLinearSystemSolver)
    {
    }

    /**
     * @brief Default destructor.
     */
    ~TrilinosChimeraBlockBuilderAndSolver() override = default;

    /**
     * Copy constructor
     */
    TrilinosChimeraBlockBuilderAndSolver(const TrilinosChimeraBlockBuilderAndSolver& rOther) = delete;

    /**
     * Assignment operator
     */
    TrilinosChimeraBlockBuilderAndSolver& operator=(const TrilinosChimeraBlockBuilderAndSolver& rOther) = delete;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void InitializeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb) override
    {
        KRATOS_TRY

        BaseType::InitializeSolutionStep(rModelPart, rA, rDx, rb);

        // Getting process info
        const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

        // Computing constraints
        const int n_constraints = static_cast<int>(rModelPart.MasterSlaveConstraints().size());
        auto constraints_begin = rModelPart.MasterSlaveConstraintsBegin();
        #pragma omp parallel for schedule(guided, 512) firstprivate(n_constraints, constraints_begin)
        for (int k = 0; k < n_constraints; ++k) {
            auto it = constraints_begin + k;
            it->InitializeSolutionStep(r_process_info); // Here each constraint constructs and stores its T and C matrices. Also its equation slave_ids.
        }

        KRATOS_CATCH("")
    }

    void FinalizeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb) override
    {
        BaseType::FinalizeSolutionStep(rModelPart, rA, rDx, rb);

        // Getting process info
        const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

        // Computing constraints
        const int n_constraints = static_cast<int>(rModelPart.MasterSlaveConstraints().size());
        const auto constraints_begin = rModelPart.MasterSlaveConstraintsBegin();
        #pragma omp parallel for schedule(guided, 512) firstprivate(n_constraints, constraints_begin)
        for (int k = 0; k < n_constraints; ++k) {
            auto it = constraints_begin + k;
            it->FinalizeSolutionStep(r_process_info);
        }
    }


    /**
     * @brief Function to perform the build the system matrix and the residual
     * vector
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param rA The LHS matrix
     * @param rb The RHS vector
     */
    void Build(typename TSchemeType::Pointer pScheme,
               ModelPart& rModelPart,
               TSystemMatrixType& rA,
               TSystemVectorType& rb) override
    {
        KRATOS_TRY
        KRATOS_ERROR_IF(!pScheme) << "No scheme provided!" << std::endl;
        BaseType::Build(pScheme, rModelPart, rA, rb);

        KRATOS_CATCH("")
    }

    /**
     * @brief Function to perform the building of the LHS
     * @details Depending on the implementation choosen the size of the matrix
     * could be equal to the total number of Dofs or to the number of
     * unrestrained dofs
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param rA The LHS matrix
     */
    void BuildLHS(typename TSchemeType::Pointer pScheme,
                  ModelPart& rModelPart,
                  TSystemMatrixType& rA) override
    {
        KRATOS_TRY
        BaseType::BuildLHS(pScheme, rModelPart, rA);
        KRATOS_CATCH("")
    }

    /**
     * @brief Build a rectangular matrix of size n*N where "n" is the number of
     * unrestrained degrees of freedom and "N" is the total number of degrees of
     * freedom involved.
     * @details This matrix is obtained by building the total matrix without the
     * lines corresponding to the fixed degrees of freedom (but keeping the
     * columns!!)
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param A The LHS matrix
     */
    void BuildLHS_CompleteOnFreeRows(typename TSchemeType::Pointer pScheme,
                                     ModelPart& rModelPart,
                                     TSystemMatrixType& A) override
    {
        KRATOS_ERROR << "Method BuildLHS_CompleteOnFreeRows not implemented in "
                        "Trilinos Builder And Solver"
                     << std::endl;
    }


    /**
     * @brief Function to perform the building and solving phase at the same time.
     * @details It is ideally the fastest and safer function to use when it is
     * possible to solve just after building
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param rA The LHS matrix
     * @param rDx The Unknowns vector
     * @param rb The RHS vector
     */
    void BuildAndSolve(typename TSchemeType::Pointer pScheme,
                       ModelPart& rModelPart,
                       TSystemMatrixType& rA,
                       TSystemVectorType& rDx,
                       TSystemVectorType& rb) override
    {
        KRATOS_TRY
        const int global_num_constraints = GetGlobalNumberOfConstraints(rModelPart);
        const Epetra_CrsGraph a_graph = rA.Graph();

        auto start_build_time = std::chrono::steady_clock::now();
        Build(pScheme, rModelPart, rA, rb);
        auto end_build_time = std::chrono::steady_clock::now();
        KRATOS_INFO_IF("TrilinosBuilderAndSolverWithConstraints",BaseType::GetEchoLevel() > 0)
                <<"Build time : "<< std::chrono::duration_cast<std::chrono::milliseconds>(end_build_time - start_build_time).count()/1000.0 <<"s"<<std::endl;
        if(global_num_constraints > 0) {
            auto start_build_ms_time = std::chrono::steady_clock::now();
            BuildMasterSlaveConstraints(rModelPart);
            auto end_build_ms_time = std::chrono::steady_clock::now();
            KRATOS_INFO_IF("TrilinosBuilderAndSolverWithConstraints",BaseType::GetEchoLevel() > 0)
                    <<"Build Master-Slave constraints time : "<< std::chrono::duration_cast<std::chrono::milliseconds>(end_build_ms_time - start_build_ms_time).count()/1000.0 <<" s"<<std::endl;

            auto start_ac_time = std::chrono::steady_clock::now();
            // ApplyConstraints(pScheme, rModelPart, rA, rb);
            auto end_ac_time = std::chrono::steady_clock::now();
            KRATOS_INFO_IF("TrilinosBuilderAndSolverWithConstraints",BaseType::GetEchoLevel() > 0)
                    <<"Apply constraints time : "<< std::chrono::duration_cast<std::chrono::milliseconds>(end_ac_time - start_ac_time).count()/1000.0 <<"s"<<std::endl;
        }

        // apply dirichlet conditions
        BaseType::ApplyDirichletConditions(pScheme, rModelPart, rA, rDx, rb);

        KRATOS_INFO_IF("TrilinosResidualBasedBlockBuilderAndSolver", BaseType::GetEchoLevel() == 3)
            << "\nBefore the solution of the system"
            << "\nSystem Matrix = " << rA << "\nunknowns vector = " << rDx
            << "\nRHS vector = " << rb << std::endl;


        auto start_solve_time = std::chrono::steady_clock::now();
        SystemSolveWithPhysics(rA, rDx, rb, rModelPart);
        auto end_solve_time = std::chrono::steady_clock::now();
        KRATOS_INFO_IF("TrilinosBuilderAndSolverWithConstraints",BaseType::GetEchoLevel() > 0)
                <<"Solve time : "<< std::chrono::duration_cast<std::chrono::milliseconds>(end_solve_time - start_solve_time).count()/1000.0 <<" s"<<std::endl;

        KRATOS_INFO_IF("TrilinosResidualBasedBlockBuilderAndSolver", BaseType::GetEchoLevel() == 3)
            << "\nAfter the solution of the system"
            << "\nSystem Matrix = " << rA << "\nUnknowns vector = " << rDx
            << "\nRHS vector = " << rb << std::endl;

        if(global_num_constraints > 0) {
            TSystemMatrixType new_a(Copy, a_graph);
            TSparseSpace::SetToZero(rA);
            TSparseSpace::Copy(new_a, rA);
        }
        KRATOS_CATCH("")
    }

    void SystemSolveWithPhysics(
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb,
        ModelPart& rModelPart
    )
    {
        const int global_num_constraints =  GetGlobalNumberOfConstraints(rModelPart);
        if(global_num_constraints > 0) {
            TSparseSpace::SetToZero(rDx);
            TSystemVectorType dx_mod(rb.Map());
            InternalSystemSolveWithPhysics(rA, dx_mod, rb, rModelPart);
            //recover solution of the original problem
            // TODO: Sparse matrix vector multiplication to get Dx
            // TSparseSpace::Mult(*mpT, dx_mod, rDx);
            rDx = dx_mod;
        } else {
            InternalSystemSolveWithPhysics(rA, rDx, rb, rModelPart);
        }
    }

    /**
     * @brief This is a call to the linear system solver
     * @param A The LHS matrix
     * @param Dx The Unknowns vector
     * @param b The RHS vector
     */
    void InternalSystemSolveWithPhysics(TSystemMatrixType &rA,
                                TSystemVectorType &rDx,
                                TSystemVectorType &rb,
                                ModelPart &rModelPart)
    {
        KRATOS_TRY

        double norm_b;
        if (TSparseSpace::Size(rb) != 0)
            norm_b = TSparseSpace::TwoNorm(rb);
        else
            norm_b = 0.00;

        if (norm_b != 0.00)
        {
            if (BaseType::mpLinearSystemSolver->AdditionalPhysicalDataIsNeeded())
                BaseType::mpLinearSystemSolver->ProvideAdditionalData(
                    rA, rDx, rb, BaseType::mDofSet, rModelPart);

            BaseType::mpLinearSystemSolver->Solve(rA, rDx, rb);
        }
        else
        {
            TSparseSpace::SetToZero(rDx);
            KRATOS_WARNING(
                "TrilinosResidualBasedBlockBuilderAndSolver")
                << "ATTENTION! setting the RHS to zero!" << std::endl;
        }

        // prints informations about the current time
        KRATOS_INFO_IF("TrilinosResidualBasedBlockBuilderAndSolver", (BaseType::GetEchoLevel() > 1))
            << *(BaseType::mpLinearSystemSolver) << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief Resizes the system matrix and the vector according to the number
     * of dos in the current rModelPart. This function also decides on the
     * sparsity pattern and the graph of the trilinos csr matrix
     * @param pScheme The integration scheme considered
     * @param rpA The LHS matrix
     * @param rpDx The Unknowns vector
     * @param rpd The RHS vector
     * @param rModelPart The model part of the problem to solve
     */
    void ResizeAndInitializeVectors(typename TSchemeType::Pointer pScheme,
                                    TSystemMatrixPointerType& rpA,
                                    TSystemVectorPointerType& rpDx,
                                    TSystemVectorPointerType& rpb,
                                    ModelPart& rModelPart) override
    {
        KRATOS_TRY
        BaseType::ResizeAndInitializeVectors(pScheme, rpA, rpDx, rpb, rModelPart);

        ConstructMasterSlaveConstraintsStructure(rModelPart);
        KRATOS_CATCH("")
    }

    /**
     * @brief This function is intended to be called at the end of the solution step to clean up memory storage not needed
     */
    void Clear() override
    {
        BaseType::Clear();
        TSparseSpace::Clear(mpL);
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

protected:
    ///@name Protected static Member Variables
    ///@{

    TSystemMatrixPointerType mpT;
    TSystemMatrixPointerType mpL;
    TSystemVectorPointerType mpConstantVector; /// This is vector containing the rigid movement of the constraint
    std::vector<IndexType> mSlaveIds;  /// The equation ids of the slaves
    std::vector<IndexType> mMasterIds; /// The equation ids of the master
    std::set<int> mInactiveSlaveEqIDs; /// The set containing the inactive slave dofs

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    void ApplyConstraints(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rb) override
    {
        const int step_num = rModelPart.GetProcessInfo().GetValue(STEP);
        const int global_num_constraints = GetGlobalNumberOfConstraints(rModelPart);

        if (global_num_constraints > 0) {
            {// To be able to clear res_b automatically.
                auto& map_b = rb.Map();
                TSystemVectorType res_b(map_b);
                const double zero = 0.0;
                TSparseSpace::TransposeMult(*mpL, rb, res_b);
                res_b.GlobalAssemble();
                // Apply diagonal values on slaves
                for (int i = 0; i < static_cast<int>(mSlaveIds.size()); ++i) {
                    const int slave_equation_id = mSlaveIds[i];
                    if (mInactiveSlaveEqIDs.find(slave_equation_id) == mInactiveSlaveEqIDs.end()) {
                        res_b.ReplaceGlobalValues(1, &slave_equation_id, &zero);
                    }
                }
                TSparseSpace::Copy(res_b, rb);
            }
            int err = 0;

            // First we do aux = T'A
            { // To delete aux_mat
                TSystemMatrixType aux_mat(Copy, mpT->RowMap(), 0);
                err = EpetraExt::MatrixMatrix::Multiply(*mpL, true, rA, false, aux_mat, false);
                KRATOS_ERROR_IF(err != 0)<<"EpetraExt MatrixMatrix multiplication(T'*A) not successful !"<<std::endl;
                aux_mat.FillComplete();
                { // To delete mod_a
                    TSystemMatrixType mod_a(Copy, aux_mat.RowMap(), 0);
                    // Now we do A = aux*T
                    // TSparseSpace::SetToZero(rA);
                    err = EpetraExt::MatrixMatrix::Multiply(aux_mat, false, *mpT, false, mod_a, false);
                    KRATOS_ERROR_IF(err != 0)<<"EpetraExt MatrixMatrix multiplication(aux*A) not successful !"<<std::endl;
                    const double inf_norm_a = rA.NormInf();
                    // Apply diagonal values on slaves
                    for (int i = 0; i < static_cast<int>(mSlaveIds.size()); ++i) {
                        const int slave_equation_id = mSlaveIds[i];
                        if (mInactiveSlaveEqIDs.find(slave_equation_id) == mInactiveSlaveEqIDs.end()) {
                            err = mod_a.ReplaceGlobalValues(slave_equation_id, 1, &inf_norm_a, &slave_equation_id);
                            if(err > 0){ // This means that the indices do not exist and we need to insert.
                                err = mod_a.InsertGlobalValues(slave_equation_id, 1, &inf_norm_a, &slave_equation_id);
                                KRATOS_ERROR_IF(err < 0)<<"Error in : InsertGlobalValues !"<<std::endl;
                            }
                        }
                    }
                    mod_a.GlobalAssemble();
                    TSparseSpace::Copy(mod_a, rA);
                }
            }
        }
    }

    virtual void BuildMasterSlaveConstraints(ModelPart& rModelPart)
    {
        KRATOS_TRY

        TSparseSpace::SetToZero(*mpT);
        TSparseSpace::SetToZero(*mpL);
        // TSparseSpace::SetToZero(mConstantVector);

        // The current process info
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Vector containing the localization in the system of the different terms
        DofsVectorType slave_dof_list, master_dof_list;

        // Contributions to the system
        Matrix transformation_matrix = LocalSystemMatrixType(0, 0);
        Vector constant_vector = LocalSystemVectorType(0);

        // Vector containing the localization in the system of the different terms
        EquationIdVectorType slave_equation_ids, master_equation_ids;

        const int number_of_constraints = static_cast<int>(rModelPart.MasterSlaveConstraints().size());
        std::set<IndexType> slave_eq_ids_set;

        // We clear the set
        mInactiveSlaveEqIDs.clear();

        for (int i_const = 0; i_const < number_of_constraints; ++i_const) {
            auto it_const = rModelPart.MasterSlaveConstraints().begin() + i_const;
            // Detect if the constraint is active or not. If the user did not make any choice the constraint
            // It is active by default
            bool constraint_is_active = it_const->IsDefined(ACTIVE) ? it_const->Is(ACTIVE) : true;
            if (constraint_is_active) {
                it_const->CalculateLocalSystem(transformation_matrix, constant_vector, r_current_process_info);
                it_const->EquationIdVector(slave_equation_ids, master_equation_ids, r_current_process_info);
                slave_eq_ids_set.insert(slave_equation_ids.begin(), slave_equation_ids.end());
                // Assemble transformation matrix
                AssembleTMatrixContribution(*mpT, transformation_matrix, slave_equation_ids, master_equation_ids);
                // Assemble the constant vector
                TSparseSpace::AssembleRHS(*mpConstantVector, constant_vector, slave_equation_ids);
            } else { // Taking into account inactive constraints
                it_const->EquationIdVector(slave_equation_ids, master_equation_ids, r_current_process_info);
                mInactiveSlaveEqIDs.insert(slave_equation_ids.begin(), slave_equation_ids.end());
            }
        }

        // All other Dofs except slaves
        double value = 1.0;
        for(auto const dof : BaseType::mDofSet){
            int eq_id = dof.EquationId();
            int my_rank = rModelPart.GetCommunicator().MyPID();
            if(my_rank == dof.GetSolutionStepValue(PARTITION_INDEX))
                if(slave_eq_ids_set.count(eq_id)==0){ // Its a master
                    mpT->ReplaceGlobalValues(1, &eq_id, 1, &eq_id, &value);
                    mpL->ReplaceGlobalValues(1, &eq_id, 1, &eq_id, &value);
                }
        }

        // All inactive slave Dofs
        // TODO: may be we should store dofs and check p index
        for(auto const inactive_id : mInactiveSlaveEqIDs){
            mpL->ReplaceGlobalValues(1, &inactive_id, 1, &inactive_id, &value);
        }

        mpT->GlobalAssemble();
        mpL->GlobalAssemble();

        TSparseSpace::WriteMatrixMarketMatrix("T_parallel.mm", *mpT, false);

        KRATOS_CATCH("")
    }

    virtual void ConstructMasterSlaveConstraintsStructure(ModelPart& rModelPart)
    {
        const int global_num_constraints = GetGlobalNumberOfConstraints(rModelPart);
        if (global_num_constraints > 0) {
            const int my_rank = rModelPart.GetCommunicator().MyPID();
            TSparseSpace::Clear(mpL);
            TSparseSpace::Clear(mpT);
            TSparseSpace::Clear(mpConstantVector);
            IndexType number_of_local_dofs = BaseType::mLastMyId - BaseType::mFirstMyId;
            std::vector<int> local_eq_ids;
            local_eq_ids.reserve(number_of_local_dofs);
            // generate map - use the "temp" array here
            // for (IndexType i = 0; i != number_of_local_dofs; i++)
            //     local_eq_ids.push_back( BaseType::mFirstMyId + i );
            for(const auto& dof : BaseType::mDofSet){
                const auto& node = rModelPart.Nodes()[dof.Id()];
                const bool is_bound = node.IsDefined(SLAVE);
                const bool is_slave =  node.IsDefined(SLAVE) ? node.Is(SLAVE) : false;
                const bool is_local = node.GetSolutionStepValue(PARTITION_INDEX) == my_rank;

                const bool is_bound_slave = is_bound && is_slave;
                const bool is_local_not_bound = is_local && !is_bound;

                if(is_bound_slave || is_local_not_bound)
                   local_eq_ids.push_back( dof.EquationId() );
            }

            std::sort(local_eq_ids.begin(), local_eq_ids.end());

            KRATOS_INFO_ALL_RANKS("T global num ")<<rModelPart.GetCommunicator().GetDataCommunicator().SumAll(local_eq_ids.size())<<std::endl;

            const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

            // Constraint initial iterator
            const auto it_const_begin = rModelPart.MasterSlaveConstraints().begin();
            std::map<IndexType, std::set<IndexType>> indices;

            Element::EquationIdVectorType slave_ids;
            Element::EquationIdVectorType master_ids;

            for (int i_const = 0; i_const < static_cast<int>(rModelPart.MasterSlaveConstraints().size()); ++i_const) {
                auto it_const = it_const_begin + i_const;

                // Detect if the constraint is active or not. If the user did not make any choice the constraint
                // It is active by default
                bool constraint_is_active = it_const->IsDefined(ACTIVE) ? it_const->Is(ACTIVE) : true;
                if(constraint_is_active) {
                    it_const->EquationIdVector(slave_ids, master_ids, r_current_process_info);
                    // Slave DoFs
                    for (auto &id_i : slave_ids) {
                        indices[id_i].insert(master_ids.begin(), master_ids.end());
                    }
                    mMasterIds.insert(mMasterIds.end(), master_ids.begin(), master_ids.end());
                }
            }

            mSlaveIds.clear();
            std::sort( mMasterIds.begin(), mMasterIds.end() );
            mMasterIds.erase( unique( mMasterIds.begin(), mMasterIds.end() ), mMasterIds.end() );
            for(const auto& slave_masters_pair : indices){
                mSlaveIds.push_back(slave_masters_pair.first);
            }
            // Diagonal values
            for(const auto& local_eq_id:local_eq_ids)
                indices[local_eq_id].insert(local_eq_id);

            // Count the row sizes
            std::size_t nnz = 10; // This is a guess. This means each slave has some number of masters
            const Epetra_Map row_map(-1, local_eq_ids.size(), local_eq_ids.data(), 0, BaseType::mrComm);

            std::set<IndexType> column_eq_ids_set;
            for(auto& slave_masters_pair : indices){
                for (auto it = slave_masters_pair.second.begin(); it != slave_masters_pair.second.end(); ++it) {
                    column_eq_ids_set.insert(*it);
                }
            }

            std::vector<int> col_eq_ids_vector(column_eq_ids_set.begin(), column_eq_ids_set.end());
            column_eq_ids_set.clear();
            const Epetra_Map col_map(-1, col_eq_ids_vector.size(), col_eq_ids_vector.data(), 0, BaseType::mrComm);
            Epetra_FECrsGraph t_graph(Copy, row_map, col_map, nnz);
            col_eq_ids_vector.clear();
            // Actually inserting indices into the graph
            int slave_eq_id = 0;
            for(auto& slave_masters_pair : indices){
                slave_eq_id = slave_masters_pair.first;
                std::vector<int> master_eq_ids(slave_masters_pair.second.begin(), slave_masters_pair.second.end());

                int ierr = t_graph.InsertGlobalIndices(1, &slave_eq_id, master_eq_ids.size(), master_eq_ids.data());
                KRATOS_ERROR_IF(ierr != 0)
                    << ": Epetra failure in Graph.InsertGlobalIndices. Error code: " << ierr
                    << std::endl;

                slave_masters_pair.second.clear(); //deallocating the memory
            }

            // The diagonal values everywhere except at the slaves
            for(const auto& eq_id : local_eq_ids)
            {
                if(indices.count(eq_id) == 0)
                {
                    int ierr = t_graph.InsertGlobalIndices(1, &eq_id, 1, &eq_id);
                    KRATOS_ERROR_IF(ierr != 0)
                        << ": Epetra failure in Graph.InsertGlobalIndices. Error code: " << ierr
                        << std::endl;
                }
            }

            int ierr = t_graph.GlobalAssemble();
            KRATOS_ERROR_IF(ierr != 0)
                << ": Epetra failure in Graph.GlobalAssemble. Error code: " << ierr
                << std::endl;

            // generate a new matrix pointer according to this graph
            TSystemMatrixPointerType p_new_t =
                TSystemMatrixPointerType(new TSystemMatrixType(Copy, t_graph));
            mpT.swap(p_new_t);

            TSystemMatrixPointerType p_new_l =
                TSystemMatrixPointerType(new TSystemMatrixType(Copy, t_graph));
            mpL.swap(p_new_l);

            TSystemVectorPointerType p_new_Cvec =
                TSystemVectorPointerType(new TSystemVectorType(row_map));
            mpConstantVector.swap(p_new_Cvec);

        }
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

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{
    int GetGlobalNumberOfConstraints(ModelPart& rModelPart){
        const auto& r_data_comm = rModelPart.GetCommunicator().GetDataCommunicator();
        const int local_num_constraints = rModelPart.MasterSlaveConstraints().size();
        const int global_num_constraints = r_data_comm.SumAll(local_num_constraints);
        return global_num_constraints;
    }


    inline static void AssembleTMatrixContribution(
        TSystemMatrixType& rMatrix,
        LocalSystemMatrixType& rContribution,
        EquationIdVectorType& rRowEquationIds,
        EquationIdVectorType& rColEquationIds
    )
    {
        //size Epetra vectors
        Epetra_IntSerialDenseVector row_indices(rRowEquationIds.size());
        Epetra_IntSerialDenseVector col_indices(rColEquationIds.size());
        Epetra_SerialDenseMatrix values(rRowEquationIds.size(), rColEquationIds.size());

        //fill epetra vectors
        int loc_i = 0;
        for (unsigned int i = 0; i < rRowEquationIds.size(); ++i)
        {
            row_indices[loc_i] = rRowEquationIds[i];

            int loc_j = 0;
            for (unsigned int j = 0; j < rColEquationIds.size(); ++j)
            {
                col_indices[loc_j] = rColEquationIds[j];
                values(loc_i, loc_j) = rContribution(i, j);
                ++loc_j;
            }
            ++loc_i;
        }

        int ierr = rMatrix.SumIntoGlobalValues(row_indices, col_indices, values);
        if(ierr != 0) KRATOS_ERROR<<"Epetra failure found"<<std::endl;
    }
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
}; /* Class TrilinosChimeraBlockBuilderAndSolver */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos.*/

#endif /* TRILINOS_CHIMERA_BLOCK_BUILDER_AND_SOLVER_WITH_CONSTRAINTS  defined */
