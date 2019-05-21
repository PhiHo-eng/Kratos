//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//


// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_conditions/grid_based_conditions/mpm_grid_line_load_condition_2d.h"
#include "utilities/math_utils.h"
#include "utilities/integration_utilities.h"

namespace Kratos
{
//******************************* CONSTRUCTOR ****************************************
//************************************************************************************

MPMGridLineLoadCondition2D::MPMGridLineLoadCondition2D( IndexType NewId, GeometryType::Pointer pGeometry )
    : MPMBaseLoadCondition( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************

MPMGridLineLoadCondition2D::MPMGridLineLoadCondition2D( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
    : MPMBaseLoadCondition( NewId, pGeometry, pProperties )
{
}

//********************************* CREATE *******************************************
//************************************************************************************

Condition::Pointer MPMGridLineLoadCondition2D::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_intrusive<MPMGridLineLoadCondition2D>(NewId, pGeom, pProperties);
}

//************************************************************************************
//************************************************************************************

Condition::Pointer MPMGridLineLoadCondition2D::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_intrusive<MPMGridLineLoadCondition2D>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

//******************************* DESTRUCTOR *****************************************
//************************************************************************************

MPMGridLineLoadCondition2D::~MPMGridLineLoadCondition2D()
{
}

//************************************************************************************
//************************************************************************************

void MPMGridLineLoadCondition2D::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo,
    bool CalculateStiffnessMatrixFlag,
    bool CalculateResidualVectorFlag
    )
{
    KRATOS_TRY;

    const GeometryType& rGeom = GetGeometry();
    const unsigned int number_of_nodes = rGeom.size();
    const unsigned int dimension = rGeom.WorkingSpaceDimension();
    const unsigned int block_size = this->GetBlockSize();

    // Resizing as needed the LHS
    unsigned int mat_size = number_of_nodes * block_size;

    if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != mat_size )
        {
            rLeftHandSideMatrix.resize( mat_size, mat_size, false );
        }

        noalias( rLeftHandSideMatrix ) = ZeroMatrix(mat_size,mat_size); //resetting LHS
    }

    //resizing as needed the RHS
    if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
    {
        if ( rRightHandSideVector.size( ) != mat_size )
        {
            rRightHandSideVector.resize( mat_size, false );
        }

        noalias( rRightHandSideVector ) = ZeroVector( mat_size ); //resetting RHS
    }

    IntegrationMethod integration_method = IntegrationUtilities::GetIntegrationMethodForExactMassMatrixEvaluation(rGeom);

    // Reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = rGeom.IntegrationPoints(integration_method);

    const GeometryType::ShapeFunctionsGradientsType& DN_De = rGeom.ShapeFunctionsLocalGradients(integration_method);

    const Matrix& Ncontainer = rGeom.ShapeFunctionsValues(integration_method);

    // Sizing work matrices
    Vector pressure_on_nodes = ZeroVector( number_of_nodes );

    // Pressure applied to the element itself
    double pressure_on_condition = 0.0;
    if( this->Has( PRESSURE ) )
    {
        pressure_on_condition += this->GetValue( PRESSURE );
    }
    if( this->Has( NEGATIVE_FACE_PRESSURE ) )
    {
        pressure_on_condition += this->GetValue( NEGATIVE_FACE_PRESSURE );
    }
    if( this->Has( POSITIVE_FACE_PRESSURE ) )
    {
        pressure_on_condition -= this->GetValue( POSITIVE_FACE_PRESSURE );
    }

    for ( unsigned int i = 0; i < pressure_on_nodes.size(); i++ )
    {
        pressure_on_nodes[i] = pressure_on_condition;
        if( rGeom[i].SolutionStepsDataHas( NEGATIVE_FACE_PRESSURE) )
        {
            pressure_on_nodes[i] += rGeom[i].FastGetSolutionStepValue( NEGATIVE_FACE_PRESSURE );
        }
        if( rGeom[i].SolutionStepsDataHas( POSITIVE_FACE_PRESSURE) )
        {
            pressure_on_nodes[i] -= rGeom[i].FastGetSolutionStepValue( POSITIVE_FACE_PRESSURE );
        }
    }

    // Vector with a loading applied to the elemnt
    array_1d<double, 3 > line_load = ZeroVector(3);
    if( this->Has( LINE_LOAD ) )
    {
        noalias(line_load) = this->GetValue( LINE_LOAD );
    }

    for ( unsigned int point_number = 0; point_number < integration_points.size(); point_number++ )
    {
        const double det_j = rGeom.DeterminantOfJacobian( integration_points[point_number] );

        const double integration_weight = GetIntegrationWeight(integration_points, point_number, det_j);

        array_1d<double, 3> normal;
        if(rGeom.WorkingSpaceDimension() == 2 )
        {
            noalias(normal) = rGeom.UnitNormal( integration_points[point_number] );
        }
        else{
            if(!Has(LOCAL_AXIS_2))
                KRATOS_ERROR << "the variable LOCAL_AXES_2 is needed to compute the normal";
            const auto& v2 = GetValue(LOCAL_AXIS_2);

            array_1d<double,3> v1 = rGeom[1].Coordinates() - rGeom[0].Coordinates();

            MathUtils<double>::CrossProduct(normal,v1,v2 );
            normal /= norm_2(normal);
        }


        // Calculating the pressure on the gauss point
        double gauss_pressure = 0.0;
        for ( unsigned int ii = 0; ii < number_of_nodes; ii++ )
        {
            gauss_pressure += Ncontainer( point_number, ii ) * pressure_on_nodes[ii];
        }

        if ( CalculateStiffnessMatrixFlag == true )
        {
            if ( std::abs(gauss_pressure) > std::numeric_limits<double>::epsilon() )
            {
                CalculateAndSubKp( rLeftHandSideMatrix, DN_De[point_number], row( Ncontainer, point_number ), gauss_pressure, integration_weight );
            }
        }
        // Adding contributions to the residual vector
        if ( CalculateResidualVectorFlag == true )
        {
            if ( std::abs(gauss_pressure) > std::numeric_limits<double>::epsilon() )
            {
                CalculateAndAddPressureForce( rRightHandSideVector, row( Ncontainer, point_number ), normal, gauss_pressure, integration_weight );
            }
        }

        array_1d<double,3> gauss_load = line_load;
        for (unsigned int ii = 0; ii < number_of_nodes; ++ii)
        {
            if( rGeom[ii].SolutionStepsDataHas( LINE_LOAD ) )
            {
                noalias(gauss_load) += ( Ncontainer( point_number, ii )) * rGeom[ii].FastGetSolutionStepValue( LINE_LOAD );
            }
        }
        for (unsigned int ii = 0; ii < number_of_nodes; ++ii)
        {
            unsigned int base = ii * block_size;

            for(unsigned int k = 0; k < dimension; ++k)
            {
                rRightHandSideVector[base + k] += integration_weight * Ncontainer( point_number, ii ) * gauss_load[k];
            }
        }
    }

    //if (this->HasRotDof()) this->CalculateAndAddWorkEquivalentNodalForcesLineLoad(gauss_load,rRightHandSideVector);

    KRATOS_CATCH( "" )
}

//***********************************************************************
//***********************************************************************

void MPMGridLineLoadCondition2D::CalculateAndSubKp(
    Matrix& K,
    const Matrix& DN_De,
    const RowMatrix& N,
    const double Pressure,
    const double IntegrationWeight
    )
{
    KRATOS_TRY

    Matrix Kij( 2, 2 );
    BoundedMatrix<double,2,2> Cross_gn( 2, 2 );

    //TODO: decide what to do with thickness
    //const double h0 = GetProperties()[THICKNESS];
    const double h0 = 1.00;
    Cross_gn( 0, 0 ) = 0.0;
    Cross_gn( 0, 1 ) = -h0;
    Cross_gn( 1, 0 ) = -h0;
    Cross_gn( 1, 1 ) = 0.0;

    for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
    {
        const unsigned int RowIndex = i * 2;

        for ( unsigned int j = 0; j < GetGeometry().size(); j++ )
        {
            const unsigned int ColIndex = j * 2;

            const double coeff = Pressure * N[i] * DN_De( j, 0 ) * IntegrationWeight;
            Kij = -coeff * Cross_gn;

            //TAKE CARE: the load correction matrix should be SUBTRACTED not added
            MathUtils<double>::SubtractMatrix( K, Kij, RowIndex, ColIndex );
        }
    }

    KRATOS_CATCH( "" )
}

//***********************************************************************
//***********************************************************************

void MPMGridLineLoadCondition2D::CalculateAndAddPressureForce(
    Vector& rRightHandSideVector,
    const RowMatrix& N,
    const array_1d<double, 3>& Normal,
    double Pressure,
    double IntegrationWeight
    )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int block_size = this->GetBlockSize();

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = block_size * i;
        const double coeff = Pressure * N[i] * IntegrationWeight;

        rRightHandSideVector[index   ]  -= coeff * Normal[0];
        rRightHandSideVector[index + 1] -= coeff * Normal[1];
    }
}

} // Namespace Kratos


