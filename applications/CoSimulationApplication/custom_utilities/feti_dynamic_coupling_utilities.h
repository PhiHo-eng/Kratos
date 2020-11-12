//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Peter Wilson
//

#if !defined(KRATOS_FETI_DYNAMIC_COUPLING_UTILITIES_H_INCLUDED)
#define  KRATOS_FETI_DYNAMIC_COUPLING_UTILITIES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "spaces/ublas_space.h"
#include "co_simulation_application_variables.h"
#include "linear_solvers/linear_solver.h"


namespace Kratos
{
    class KRATOS_API(CO_SIMULATION_APPLICATION) FetiDynamicCouplingUtilities
    {
    public:
        typedef std::size_t SizeType;
        typedef std::size_t IndexType;

        typedef Node<3> NodeType;
        typedef typename NodeType::Pointer NodePointerType;
        typedef Geometry<NodeType> GeometryType;
        typedef typename GeometryType::Pointer GeometryPointerType;

        typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SparseSpaceType;
        typedef typename SparseSpaceType::MatrixType SystemMatrixType;

        typedef UblasSpace<double, Matrix, Vector> DenseSpaceType;

        typedef Matrix DenseMappingMatrixType;

        typedef typename SparseSpaceType::MatrixType MappingMatrixType;

        typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
        typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;
        typedef Kratos::shared_ptr<LinearSolverType> LinearSolverSharedPointerType;

        /// The definition of the numerical limit
        static constexpr double numerical_limit = std::numeric_limits<double>::epsilon();

        enum class SolverIndex { Origin, Destination };
        enum class EquilibriumVariable { Displacement, Velocity, Acceleration};

        FetiDynamicCouplingUtilities(ModelPart & rInterfaceOrigin, ModelPart & rInterFaceDestination,
            const Parameters JsonParameters);

        void SetOriginAndDestinationDomainsWithInterfaceModelParts(ModelPart& rInterfaceOrigin,
            ModelPart& rInterFaceDestination);

        void SetEffectiveStiffnessMatrixImplicit(SystemMatrixType& rK, const IndexType SolverIndex);

        void SetEffectiveStiffnessMatrixExplicit(const IndexType SolverIndex)
        {
            if (SolverIndex == 0) mSubTimestepIndex = 1;
        };

        void SetMappingMatrix(CompressedMatrix& rMappingMatrix)
        {
            mpMappingMatrix = &rMappingMatrix;
        };

        void SetLinearSolver(LinearSolverSharedPointerType pSolver)
        {
            mpSolver = pSolver;
        }

        void SetOriginInitialKinematics();

        void EquilibrateDomains();

    private:
        ModelPart& mrOriginInterfaceModelPart;
        ModelPart& mrDestinationInterfaceModelPart;

        ModelPart* mpOriginDomain = nullptr;
        ModelPart* mpDestinationDomain = nullptr;

        SystemMatrixType* mpKOrigin = nullptr;
        SystemMatrixType* mpKDestination = nullptr;

        CompressedMatrix* mpMappingMatrix = nullptr;
        CompressedMatrix* mpMappingMatrixForce = nullptr;

        // Origin quantities
        Vector mInitialOriginInterfaceKinematics;
        Vector mFinalOriginInterfaceKinematics;
        CompressedMatrix mProjectorOrigin;
        CompressedMatrix mUnitResponseOrigin;

        // Quantities to store for a linear system
        CompressedMatrix mCondensationMatrix;
        CompressedMatrix mUnitResponseDestination;
        CompressedMatrix mProjectorDestination;
        bool mIsLinearSetupComplete = false;

        EquilibriumVariable mEquilibriumVariable;

        LinearSolverSharedPointerType mpSolver = nullptr;

        bool mIsImplicitOrigin;
        bool mIsImplicitDestination;
        const Parameters mParameters;
        bool mIsLinear = false;

        IndexType mSubTimestepIndex = 1;
        IndexType mTimestepRatio;

        const bool mIsCheckEquilibrium = true; // normally true

        void CalculateUnbalancedInterfaceFreeKinematics(Vector& rUnbalancedKinematics, const bool IsEquilibriumCheck = false);

        void GetInterfaceQuantity(ModelPart& rInterface, const Variable< array_1d<double, 3> >& rVariable,
            Vector& rContainer, const SizeType nDOFs);

        void GetInterfaceQuantity(ModelPart& rInterface, const Variable<double>& rVariable,
            Vector& rContainer, const SizeType nDOFs);

        void GetExpandedMappingMatrix(CompressedMatrix& rExpandedMappingMat, const SizeType nDOFs);

        void ComposeProjector(CompressedMatrix& rProjector, const SolverIndex solverIndex);

        void DetermineDomainUnitAccelerationResponse(SystemMatrixType* pK,
            const CompressedMatrix& rProjector, CompressedMatrix& rUnitResponse, const SolverIndex solverIndex);

        void DetermineDomainUnitAccelerationResponseExplicit(CompressedMatrix& rUnitResponse,
            const CompressedMatrix& rProjector, ModelPart& rDomain, const SolverIndex solverIndex);

        void DetermineDomainUnitAccelerationResponseImplicit(CompressedMatrix& rUnitResponse,
            const CompressedMatrix& rProjector, SystemMatrixType* pK, const SolverIndex solverIndex);

        void CalculateCondensationMatrix(CompressedMatrix& rCondensationMatrix,
            const CompressedMatrix& rOriginUnitResponse, const CompressedMatrix& rDestinationUnitResponse,
            const CompressedMatrix& rOriginProjector, const CompressedMatrix& rDestinationProjector);

        void DetermineLagrangianMultipliers(Vector& rLagrangeVec,
            CompressedMatrix& rCondensationMatrix, Vector& rUnbalancedKinematics);

        void ApplyCorrectionQuantities(const Vector& rLagrangeVec,
            const CompressedMatrix& rUnitResponse, const SolverIndex solverIndex);

        void AddCorrectionToDomain(ModelPart* pDomain,
            const Variable< array_1d<double, 3> >& rVariable,
            const Vector& rCorrection, const bool IsImplicit);

        void WriteLagrangeMultiplierResults(const Vector& rLagrange);

        void ApplyMappingMatrixToProjector(CompressedMatrix& rProjector, const SizeType DOFs);

        void PrintInterfaceKinematics(const Variable< array_1d<double, 3> >& rVariable, const SolverIndex solverIndex);

        Variable< array_1d<double, 3> >& GetEquilibriumVariable();

    };  // namespace FetiDynamicCouplingUtilities.

}  // namespace Kratos.

#endif // KRATOS_FETI_DYNAMIC_COUPLING_UTILITIES_H_INCLUDED  defined
