//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Veronika Singer
//


// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_conditions/particle_based_conditions/mpm_particle_lagrange_dirichlet_condition.h"
#include "includes/kratos_flags.h"
#include "utilities/math_utils.h"
#include "custom_utilities/particle_mechanics_math_utilities.h"
#include "includes/checks.h"

namespace Kratos
{
//******************************* CONSTRUCTOR ****************************************
//************************************************************************************

MPMParticleLagrangeDirichletCondition::MPMParticleLagrangeDirichletCondition( IndexType NewId, GeometryType::Pointer pGeometry )
    : MPMParticleBaseDirichletCondition( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************

MPMParticleLagrangeDirichletCondition::MPMParticleLagrangeDirichletCondition( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
    : MPMParticleBaseDirichletCondition( NewId, pGeometry, pProperties )
{
}

//********************************* CREATE *******************************************
//************************************************************************************

Condition::Pointer MPMParticleLagrangeDirichletCondition::Create(IndexType NewId,GeometryType::Pointer pGeometry,PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<MPMParticleLagrangeDirichletCondition>(NewId, pGeometry, pProperties);
}

//************************************************************************************
//************************************************************************************

Condition::Pointer MPMParticleLagrangeDirichletCondition::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<MPMParticleLagrangeDirichletCondition>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

//******************************* DESTRUCTOR *****************************************
//************************************************************************************

MPMParticleLagrangeDirichletCondition::~MPMParticleLagrangeDirichletCondition()
{
}

//************************************************************************************
//************************************************************************************

void MPMParticleLagrangeDirichletCondition::InitializeSolutionStep( const ProcessInfo& rCurrentProcessInfo )
{
    MPMParticleBaseDirichletCondition::InitializeSolutionStep( rCurrentProcessInfo );

    GeometryType& r_geometry = GetGeometry();

    const unsigned int number_of_nodes = r_geometry.PointsNumber();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    auto pBoundaryParticle = GetValue(MPC_LAGRANGE_NODE);

    array_1d<double, 3 > & r_lagrange_multiplier  = pBoundaryParticle->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER);

    for ( unsigned int j = 0; j < dimension; j++ )
    {
        r_lagrange_multiplier[0] *= 0.0;
    }


}

//************************************************************************************
//************************************************************************************

void MPMParticleLagrangeDirichletCondition::CalculateAll(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    bool CalculateStiffnessMatrixFlag,
    bool CalculateResidualVectorFlag
    )
{
    KRATOS_TRY

    GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int block_size = this->GetBlockSize();

    // Resizing as needed the LHS
    const unsigned int matrix_size = number_of_nodes * dimension + dimension;

    if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != matrix_size )
        {
            rLeftHandSideMatrix.resize( matrix_size, matrix_size, false );
        }

        noalias( rLeftHandSideMatrix ) = ZeroMatrix(matrix_size,matrix_size); //resetting LHS
    }

    // Resizing as needed the RHS
    if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
    {
        if ( rRightHandSideVector.size( ) != matrix_size )
        {
            rRightHandSideVector.resize( matrix_size, false );
        }

        noalias( rRightHandSideVector ) = ZeroVector( matrix_size ); //resetting RHS
    }
    
    // Prepare variables
    GeneralVariables Variables;

    // Calculating shape function
    MPMShapeFunctionPointValues(Variables.N);
    Variables.CurrentDisp = this->CalculateCurrentDisp(Variables.CurrentDisp, rCurrentProcessInfo);


    bool apply_constraints = true;

    // Check nodal mass to avoid singularity; only apply condition if particles are inside
    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        if (r_geometry[i].FastGetSolutionStepValue(NODAL_MASS) <= std::numeric_limits<double>::epsilon() ) 
        {
            apply_constraints=false;      
        }
    }  

    if (Is(CONTACT))
    {
        // NOTE: the unit_normal_vector is assumed always pointing outside the boundary
        array_1d<double, 3 > field_displacement = ZeroVector(3);
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            for ( unsigned int j = 0; j < dimension; j++)
            {
                field_displacement[j] += Variables.N[i] * Variables.CurrentDisp(i,j);
            }
        }

        const double penetration = MathUtils<double>::Dot((field_displacement - m_imposed_displacement), m_unit_normal);
        // If penetrates, apply constraint, otherwise no
        if (penetration >= 0.0)
        {
            apply_constraints = false;

        }

    }
    

    if (apply_constraints)
    {
        Matrix lagrange_matrix = ZeroMatrix(matrix_size, matrix_size);

        // Loop over shape functions of displacements
        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            const unsigned int ibase = dimension * number_of_nodes;
            for (unsigned int k = 0; k < dimension; k++)
            {
                lagrange_matrix(i* dimension+k, ibase+k) = Variables.N[i];
                lagrange_matrix(ibase+k, i*dimension + k) = Variables.N[i];
    }
        }

        // Calculate LHS Matrix and RHS Vector
        if ( CalculateStiffnessMatrixFlag == true )
        {
            rLeftHandSideMatrix = lagrange_matrix * this->GetIntegrationWeight();
        }

        if ( CalculateResidualVectorFlag == true )
        {
            //last rows of RHS
            Vector gap_function = ZeroVector(matrix_size);
            Vector right_hand_side = ZeroVector(matrix_size);
            for (unsigned int i = 0; i < number_of_nodes; i++)
            {
                const array_1d<double, 3>& r_displacement = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT);
                const int index = dimension * i;

                for (unsigned int j = 0; j < dimension; j++)
                {
                    gap_function[index+j] = r_displacement[j] ;
                }


            }
            right_hand_side = prod(lagrange_matrix, gap_function);

            //first rows of RHS
            gap_function = ZeroVector(matrix_size);
            auto pBoundaryParticle = GetValue(MPC_LAGRANGE_NODE);
            const array_1d<double, 3>& r_lagrange_multiplier = pBoundaryParticle->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER);

            for (unsigned int j = 0; j < dimension; j++){
                gap_function[dimension * number_of_nodes+j] = r_lagrange_multiplier[j];
            }
            
            right_hand_side += prod(lagrange_matrix, gap_function);

            //add imposed displacement
            for (unsigned int j = 0; j < dimension; j++){
                right_hand_side[dimension * number_of_nodes+j] -= m_imposed_displacement[j];
            }

            

            right_hand_side *= this->GetIntegrationWeight();
            noalias(rRightHandSideVector) = -right_hand_side;
            
        }
        }
    else{

        KRATOS_WATCH("NO CONDITIONS APPLIED")

    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void MPMParticleLagrangeDirichletCondition::FinalizeSolutionStep( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    MPMParticleBaseDirichletCondition::FinalizeSolutionStep(rCurrentProcessInfo);
    auto pBoundaryParticle = GetValue(MPC_LAGRANGE_NODE);
    const array_1d<double, 3>& r_lagrange_multiplier = pBoundaryParticle->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER);
    m_contact_force = r_lagrange_multiplier;

    KRATOS_CATCH( "" )
}

void MPMParticleLagrangeDirichletCondition::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY

    const GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size() ;
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    if (rResult.size() != dimension * number_of_nodes + dimension)
    {
        rResult.resize(dimension * number_of_nodes + dimension ,false);
    }


    for (unsigned int i = 0; i < number_of_nodes; ++i)
    {
        unsigned int index = i * dimension;
        rResult[index    ] = r_geometry[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[index + 1] = r_geometry[i].GetDof(DISPLACEMENT_Y).EquationId();
        if(dimension == 3)
            rResult[index + 2] = r_geometry[i].GetDof(DISPLACEMENT_Z).EquationId();

    }

    unsigned int index = number_of_nodes * dimension;
    auto pBoundaryParticle = GetValue(MPC_LAGRANGE_NODE);

    rResult[index    ] = pBoundaryParticle->GetDof(VECTOR_LAGRANGE_MULTIPLIER_X).EquationId();
    rResult[index + 1] = pBoundaryParticle->GetDof(VECTOR_LAGRANGE_MULTIPLIER_Y).EquationId();
    if(dimension == 3)
        rResult[index + 2] = pBoundaryParticle->GetDof(VECTOR_LAGRANGE_MULTIPLIER_Z).EquationId();


    KRATOS_CATCH("")
}

void MPMParticleLagrangeDirichletCondition::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY

    const GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    const unsigned int dimension =  r_geometry.WorkingSpaceDimension();
    rElementalDofList.resize(0);
    rElementalDofList.reserve(dimension * number_of_nodes + dimension);

    for (unsigned int i = 0; i < number_of_nodes; ++i)
    {
        rElementalDofList.push_back( r_geometry[i].pGetDof(DISPLACEMENT_X));
        rElementalDofList.push_back( r_geometry[i].pGetDof(DISPLACEMENT_Y));
        if(dimension == 3)
            rElementalDofList.push_back( r_geometry[i].pGetDof(DISPLACEMENT_Z));

    }
    auto pBoundaryParticle = GetValue(MPC_LAGRANGE_NODE);

    rElementalDofList.push_back(pBoundaryParticle->pGetDof(VECTOR_LAGRANGE_MULTIPLIER_X));
    rElementalDofList.push_back(pBoundaryParticle->pGetDof(VECTOR_LAGRANGE_MULTIPLIER_Y));
    if(dimension == 3)
        rElementalDofList.push_back(pBoundaryParticle->pGetDof(VECTOR_LAGRANGE_MULTIPLIER_Z));



    KRATOS_CATCH("")
}

void MPMParticleLagrangeDirichletCondition::GetValuesVector(
    Vector& rValues,
    int Step
    ) const
{
    const GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    const unsigned int mat_size = number_of_nodes * dimension + dimension;

    if (rValues.size() != mat_size)
    {
        rValues.resize(mat_size, false);
    }
    rValues = ZeroVector(mat_size);

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        const array_1d<double, 3 > & r_displacement = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT, Step);

        unsigned int index = i * dimension;

        for(unsigned int k = 0; k < dimension; ++k)
        {
            rValues[index + k] = r_displacement[k];
        }
    }
    auto pBoundaryParticle = GetValue(MPC_LAGRANGE_NODE);
    const array_1d<double, 3 > & r_lagrange_multiplier = pBoundaryParticle->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER, Step);
    const unsigned int lagrange_index = number_of_nodes * dimension;
    for(unsigned int k = 0; k < dimension; ++k)
        {
            rValues[lagrange_index + k] = r_lagrange_multiplier[k];
        }
}

void MPMParticleLagrangeDirichletCondition::GetFirstDerivativesVector(
    Vector& rValues,
    int Step
    ) const
{
    const GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    const unsigned int mat_size = number_of_nodes * dimension + dimension;

    if (rValues.size() != mat_size)
    {
        rValues.resize(mat_size, false);
    }
    rValues = ZeroVector(mat_size);

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        const array_1d<double, 3 > & r_velocity = r_geometry[i].FastGetSolutionStepValue(VELOCITY, Step);
        const unsigned int index = i * dimension;
        for(unsigned int k = 0; k < dimension; ++k)
        {
            rValues[index + k] = r_velocity[k];
        }
    }
}

void MPMParticleLagrangeDirichletCondition::GetSecondDerivativesVector(
    Vector& rValues,
    int Step
    ) const
{
    const GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    const unsigned int mat_size = number_of_nodes * dimension + dimension;

    if (rValues.size() != mat_size)
    {
        rValues.resize(mat_size, false);
    }
    rValues = ZeroVector(mat_size);

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        const array_1d<double, 3 > & r_acceleration = r_geometry[i].FastGetSolutionStepValue(ACCELERATION, Step);
        const unsigned int index = i * dimension;
        for(unsigned int k = 0; k < dimension; ++k)
        {
            rValues[index + k] = r_acceleration[k];
        }
    }
}

void MPMParticleLagrangeDirichletCondition::CalculateOnIntegrationPoints(const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1)
        rValues.resize(1);

    else {
        MPMParticleBaseDirichletCondition::CalculateOnIntegrationPoints(
            rVariable, rValues, rCurrentProcessInfo);
    }
}

void MPMParticleLagrangeDirichletCondition::CalculateOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
    std::vector<array_1d<double, 3 > >& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1)
        rValues.resize(1);

    if (rVariable == MPC_IMPOSED_DISPLACEMENT) {
        rValues[0] = m_imposed_displacement;
    }
    else if (rVariable == MPC_NORMAL) {
        rValues[0] = m_unit_normal;
    }
    else if (rVariable == MPC_CONTACT_FORCE) {
        rValues[0] = m_contact_force;
    }
    else {
        MPMParticleBaseDirichletCondition::CalculateOnIntegrationPoints(
            rVariable, rValues, rCurrentProcessInfo);
    }
}

void MPMParticleLagrangeDirichletCondition::SetValuesOnIntegrationPoints(const Variable<double>& rVariable,
    const std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR_IF(rValues.size() > 1)
        << "Only 1 value per integration point allowed! Passed values vector size: "
        << rValues.size() << std::endl;

    else {
        MPMParticleBaseDirichletCondition::SetValuesOnIntegrationPoints(
            rVariable, rValues, rCurrentProcessInfo);
    }
}

void MPMParticleLagrangeDirichletCondition::SetValuesOnIntegrationPoints(
    const Variable<array_1d<double, 3 > >& rVariable,
    const std::vector<array_1d<double, 3 > >& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR_IF(rValues.size() > 1)
        << "Only 1 value per integration point allowed! Passed values vector size: "
        << rValues.size() << std::endl;

    if (rVariable == MPC_IMPOSED_DISPLACEMENT) {
        m_imposed_displacement = rValues[0];
    }
    else if (rVariable == MPC_NORMAL) {
        m_unit_normal = rValues[0];
        ParticleMechanicsMathUtilities<double>::Normalize(m_unit_normal);
    }
    else if (rVariable == MPC_CONTACT_FORCE) {
        m_contact_force = rValues[0];
    }
    else {
        MPMParticleBaseDirichletCondition::SetValuesOnIntegrationPoints(
            rVariable, rValues, rCurrentProcessInfo);
    }
}

int MPMParticleLagrangeDirichletCondition::Check( const ProcessInfo& rCurrentProcessInfo ) const
{
    MPMParticleBaseDirichletCondition::Check(rCurrentProcessInfo);

    // Verify that the dofs exist
    for (const auto& r_node : this->GetGeometry().Points())
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(NORMAL,r_node)

    return 0;
}

} // Namespace Kratos

