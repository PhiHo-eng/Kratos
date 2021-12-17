//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Reza Najian Asl
//

// System includes


// External includes


// Project includes
#include "custom_elements/helmholtz_vec_element.h"
#include "shape_optimization_application_variables.h"
#include "includes/checks.h"
#include "includes/define.h"
#include "utilities/math_utils.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************
HelmholtzVecElement::HelmholtzVecElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!

}

//************************************************************************************
//************************************************************************************
HelmholtzVecElement::HelmholtzVecElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
}

Element::Pointer HelmholtzVecElement::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<HelmholtzVecElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

Element::Pointer HelmholtzVecElement::Create(IndexType NewId, GeometryType::Pointer pGeom,  PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<HelmholtzVecElement>(NewId, pGeom, pProperties);
}

HelmholtzVecElement::~HelmholtzVecElement()
{
}

//************************************************************************************
//************************************************************************************
void HelmholtzVecElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                            VectorType& rRightHandSideVector,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    auto& r_geometry = this->GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();

    // Resizing as needed the LHS
    const SizeType mat_size = number_of_nodes * dimension;

    if ( rLeftHandSideMatrix.size1() != mat_size )
        rLeftHandSideMatrix.resize( mat_size, mat_size, false );

    noalias( rLeftHandSideMatrix ) = ZeroMatrix( mat_size, mat_size ); //resetting LHS

    // Resizing as needed the RHS 
    if ( rRightHandSideVector.size() != mat_size )
        rRightHandSideVector.resize( mat_size, false );

    rRightHandSideVector = ZeroVector( mat_size ); //resetting RHS

    MatrixType M;
    CalculateBulkMassMatrix(M,rCurrentProcessInfo);
    MatrixType A;
    CalculateBulkStiffnessMatrix(A,rCurrentProcessInfo);

    MatrixType K = A + M;

    const unsigned int number_of_points = r_geometry.size();
    Vector nodal_vals(number_of_points*3);
    for(unsigned int node_element = 0; node_element<number_of_points; node_element++)
    {
        const VectorType &source = r_geometry[node_element].FastGetSolutionStepValue(HELMHOLTZ_SOURCE);
        auto node_weight = r_geometry[node_element].GetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS);
        if(rCurrentProcessInfo[COMPUTE_CONTROL_POINTS])
            node_weight = 1.0;
        nodal_vals[3 * node_element + 0] = source[0]/node_weight;
        nodal_vals[3 * node_element + 1] = source[1]/node_weight;
        nodal_vals[3 * node_element + 2] = source[2]/node_weight;
    }

    if(rCurrentProcessInfo[COMPUTE_CONTROL_POINTS]){
        noalias(rLeftHandSideMatrix) += M;
        noalias(rRightHandSideVector) += prod(K,nodal_vals);
        //apply drichlet BC
        Vector temp(number_of_points*3);
        for (SizeType iNode = 0; iNode < number_of_points; ++iNode) {
            const VectorType &vars = r_geometry[iNode].FastGetSolutionStepValue(HELMHOLTZ_VARS,0);
            temp[3*iNode] = vars[0];
            temp[3*iNode+1] = vars[1];
            temp[3*iNode+2] = vars[2];
        }
        noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp);     
    }            
    else{
        noalias(rLeftHandSideMatrix) += K;
        noalias(rRightHandSideVector) += nodal_vals;
    } 

    KRATOS_CATCH("")
}

//******************************************************************************
//******************************************************************************
void HelmholtzVecElement::GetValuesVector(VectorType &rValues,
                                            int Step) const {
  const GeometryType &rgeom = this->GetGeometry();
  const SizeType num_nodes = rgeom.PointsNumber();
  const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
  const unsigned int local_size = num_nodes * dimension;

  if (rValues.size() != local_size)
    rValues.resize(local_size, false);

  if (dimension == 2) {
    SizeType index = 0;
    for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
      rValues[index++] =
          rgeom[i_node].FastGetSolutionStepValue(HELMHOLTZ_VARS_X, Step);
      rValues[index++] =
          rgeom[i_node].FastGetSolutionStepValue(HELMHOLTZ_VARS_Y, Step);
    }
  } else if (dimension == 3) {
    SizeType index = 0;
    for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
      rValues[index++] =
          rgeom[i_node].FastGetSolutionStepValue(HELMHOLTZ_VARS_X, Step);
      rValues[index++] =
          rgeom[i_node].FastGetSolutionStepValue(HELMHOLTZ_VARS_Y, Step);
      rValues[index++] =
          rgeom[i_node].FastGetSolutionStepValue(HELMHOLTZ_VARS_Z, Step);
    }
  }
}

/***********************************************************************************/
/***********************************************************************************/

void HelmholtzVecElement::Calculate(const Variable<Matrix>& rVariable, Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == HELMHOLTZ_MASS_MATRIX)
        CalculateBulkMassMatrix(rOutput,rCurrentProcessInfo);

}
//************************************************************************************
//************************************************************************************
void HelmholtzVecElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    VectorType temp(0);
    CalculateLocalSystem(rLeftHandSideMatrix, temp, rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************
void HelmholtzVecElement::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType temp(0,0);
    CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************
void HelmholtzVecElement::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    if (rResult.size() != dimension * number_of_nodes)
        rResult.resize(dimension * number_of_nodes,false);

    const SizeType pos = this->GetGeometry()[0].GetDofPosition(HELMHOLTZ_VARS_X);

    if(dimension == 2) {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const SizeType index = i * 2;
            rResult[index] = GetGeometry()[i].GetDof(HELMHOLTZ_VARS_X,pos).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(HELMHOLTZ_VARS_Y,pos+1).EquationId();
        }
    } else {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const SizeType index = i * 3;
            rResult[index] = GetGeometry()[i].GetDof(HELMHOLTZ_VARS_X,pos).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(HELMHOLTZ_VARS_Y,pos+1).EquationId();
            rResult[index + 2] = GetGeometry()[i].GetDof(HELMHOLTZ_VARS_Z,pos+2).EquationId();
        }
    }

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void HelmholtzVecElement::GetDofList(DofsVectorType& rElementalDofList,const ProcessInfo& rCurrentProcessInfo) const
{

    KRATOS_TRY;

    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    rElementalDofList.resize(0);
    rElementalDofList.reserve(dimension*number_of_nodes);

    if(dimension == 2) {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(HELMHOLTZ_VARS_X));
            rElementalDofList.push_back( GetGeometry()[i].pGetDof(HELMHOLTZ_VARS_Y));
        }
    } else {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(HELMHOLTZ_VARS_X));
            rElementalDofList.push_back( GetGeometry()[i].pGetDof(HELMHOLTZ_VARS_Y));
            rElementalDofList.push_back( GetGeometry()[i].pGetDof(HELMHOLTZ_VARS_Z));
        }
    }

    KRATOS_CATCH("")

}

//************************************************************************************
//************************************************************************************
int HelmholtzVecElement::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    int check = Element::Check(rCurrentProcessInfo);

    const auto& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for ( IndexType i = 0; i < number_of_nodes; i++ ) {
        const NodeType &rnode = r_geometry[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(HELMHOLTZ_VARS,rnode)

        KRATOS_CHECK_DOF_IN_NODE(HELMHOLTZ_VARS_X, rnode)
        KRATOS_CHECK_DOF_IN_NODE(HELMHOLTZ_VARS_Y, rnode)
        KRATOS_CHECK_DOF_IN_NODE(HELMHOLTZ_VARS_Z, rnode)
    }

    return check;

    KRATOS_CATCH( "" );
}
/***********************************************************************************/
/***********************************************************************************/

void HelmholtzVecElement::CalculateBulkMassMatrix(
    MatrixType& rMassMatrix,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    const auto& r_geom = GetGeometry();
    SizeType dimension = r_geom.WorkingSpaceDimension();
    SizeType number_of_nodes = r_geom.size();
    SizeType mat_size = dimension * number_of_nodes;

    // Clear matrix
    if (rMassMatrix.size1() != mat_size || rMassMatrix.size2() != mat_size)
        rMassMatrix.resize( mat_size, mat_size, false );
    rMassMatrix = ZeroMatrix( mat_size, mat_size );


    Matrix J0(dimension, dimension);

    IntegrationMethod integration_method = IntegrationUtilities::GetIntegrationMethodForExactMassMatrixEvaluation(r_geom);
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geom.IntegrationPoints( integration_method );
    const Matrix& Ncontainer = r_geom.ShapeFunctionsValues(integration_method);

    for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
        GeometryUtils::JacobianOnInitialConfiguration(
            r_geom, integration_points[point_number], J0);
        const double detJ0 = MathUtils<double>::Det(J0);
        const double integration_weight = integration_points[point_number].Weight() * detJ0;
        const Vector& rN = row(Ncontainer,point_number);

        for ( IndexType i = 0; i < number_of_nodes; ++i ) {
            const SizeType index_i = i * dimension;

            for ( IndexType j = 0; j < number_of_nodes; ++j ) {
                const SizeType index_j = j * dimension;
                const double NiNj_weight = rN[i] * rN[j] * integration_weight;

                for ( IndexType k = 0; k < dimension; ++k )
                    rMassMatrix( index_i + k, index_j + k ) += NiNj_weight;
            }
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void HelmholtzVecElement::CalculateBulkStiffnessMatrix(
    MatrixType& rStiffnessMatrix,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    const auto& r_prop = GetProperties();

    // Checking radius
    KRATOS_ERROR_IF_NOT(r_prop.Has(HELMHOLTZ_RADIUS)) << "HELMHOLTZ_RADIUS has to be provided for the calculations of the HelmholtzVecElement!" << std::endl;

    const auto& r_geom = GetGeometry();
    SizeType dimension = r_geom.WorkingSpaceDimension();
    SizeType number_of_nodes = r_geom.size();
    SizeType mat_size = dimension * number_of_nodes;

    // Clear matrix
    if (rStiffnessMatrix.size1() != mat_size || rStiffnessMatrix.size2() != mat_size)
        rStiffnessMatrix.resize( mat_size, mat_size, false );
    rStiffnessMatrix = ZeroMatrix( mat_size, mat_size );

    //reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geom.IntegrationPoints(r_geom.GetDefaultIntegrationMethod());

    Element::GeometryType::JacobiansType J0;
    Matrix DN_DX(number_of_nodes,dimension);
    Matrix InvJ0(dimension,dimension);
    r_geom.Jacobian(J0,r_geom.GetDefaultIntegrationMethod());
    double DetJ0;

    for(std::size_t i_point = 0; i_point<integration_points.size(); ++i_point)
    {
        //calculating inverse jacobian and jacobian determinant
        MathUtils<double>::InvertMatrix(J0[i_point],InvJ0,DetJ0);

        MatrixType B = CalculateBMatrix(dimension, i_point);

        const double r_helmholtz = r_prop[HELMHOLTZ_RADIUS];
        MatrixType constitutive_matrix = SetAndModifyConstitutiveLaw(dimension, i_point);
        const double IntToReferenceWeight = integration_points[i_point].Weight() * DetJ0;

        noalias(rStiffnessMatrix) += r_helmholtz * r_helmholtz * prod(trans(B), IntToReferenceWeight * Matrix(prod(constitutive_matrix, B)));
        
    }

    KRATOS_CATCH("");
}

//******************************************************************************
//******************************************************************************
HelmholtzVecElement::MatrixType
HelmholtzVecElement::SetAndModifyConstitutiveLaw(
    const int Dimension, const int PointNumber) const {
  KRATOS_TRY;

  GeometryType::JacobiansType J0;
  GeometryType::JacobiansType invJ0;
  VectorType detJ0;
  const GeometryType &rgeom = this->GetGeometry();
  const IntegrationMethod this_integration_method =
      rgeom.GetDefaultIntegrationMethod();


  const auto& r_integration_points = rgeom.IntegrationPoints(this_integration_method);
  if (invJ0.size() != r_integration_points.size()) {
    invJ0.resize(r_integration_points.size());
  }
  if (detJ0.size() != r_integration_points.size()) {
    detJ0.resize(r_integration_points.size());
  }

  J0 = GetGeometry().Jacobian(J0, this_integration_method);

  MathUtils<double>::InvertMatrix(J0[PointNumber], invJ0[PointNumber],
                                  detJ0[PointNumber]);

  // Stiffening of elements using Jacobian determinants and exponent between
  // 0.0 and 2.0
  const double factor =
      100;               // Factor influences how far the HELMHOLTZ_VARS spreads
                         // into the fluid mesh
  const double xi = 1.5; // 1.5 Exponent influences stiffening of smaller
                         // elements; 0 = no stiffening
  const double quotient = factor / detJ0[PointNumber];
  const double weighting_factor = detJ0[PointNumber] * std::pow(quotient, xi);
  const double poisson_coefficient = this->pGetProperties()->Has(HELMHOLTZ_POISSON_RATIO)
    ? this->pGetProperties()->GetValue(HELMHOLTZ_POISSON_RATIO) : 0.3;

  // The ratio between lambda and mu affects relative stiffening against
  // volume or shape change.
  const double lambda =
      weighting_factor * poisson_coefficient /
      ((1 + poisson_coefficient) * (1 - 2 * poisson_coefficient));
  const double mu = weighting_factor / (2 * (1 + poisson_coefficient));

  MatrixType constitutive_matrix;

  // stress = lambda*tr(strain tensor)*I + 2*mu*(strain tensor).
  if (Dimension == 2) {
    constitutive_matrix = ZeroMatrix(3, 3);
    constitutive_matrix(0, 0) = lambda + 2 * mu;
    constitutive_matrix(1, 1) = constitutive_matrix(0, 0);
    constitutive_matrix(2, 2) = mu;
    constitutive_matrix(0, 1) = lambda;
    constitutive_matrix(1, 0) = lambda;
  }

  else if (Dimension == 3) {
    constitutive_matrix = ZeroMatrix(6, 6);
    constitutive_matrix(0, 0) = lambda + 2 * mu;
    constitutive_matrix(1, 1) = constitutive_matrix(0, 0);
    constitutive_matrix(2, 2) = constitutive_matrix(0, 0);
    constitutive_matrix(3, 3) = mu;
    constitutive_matrix(4, 4) = mu;
    constitutive_matrix(5, 5) = mu;
    constitutive_matrix(0, 1) = lambda;
    constitutive_matrix(1, 0) = lambda;
    constitutive_matrix(0, 2) = lambda;
    constitutive_matrix(2, 0) = lambda;
    constitutive_matrix(1, 2) = lambda;
    constitutive_matrix(2, 1) = lambda;
  }

  return constitutive_matrix;

  KRATOS_CATCH("");
}

//******************************************************************************
//******************************************************************************
HelmholtzVecElement::MatrixType
HelmholtzVecElement::CalculateBMatrix(const int Dimension,
                                              const int PointNumber) const {
  KRATOS_TRY;

  const GeometryType &rgeom = this->GetGeometry();
  const IntegrationMethod this_integration_method =
      rgeom.GetDefaultIntegrationMethod();
  GeometryType::ShapeFunctionsGradientsType DN_De =
      rgeom.ShapeFunctionsLocalGradients(this_integration_method);
  GeometryType::JacobiansType J0;
  GeometryType::JacobiansType invJ0;
  VectorType detJ0;

  const auto& r_integration_points = rgeom.IntegrationPoints(this_integration_method);
  if (invJ0.size() != r_integration_points.size()) {
    invJ0.resize(r_integration_points.size());
  }
  if (detJ0.size() != r_integration_points.size()) {
    detJ0.resize(r_integration_points.size());
  }

  J0 = GetGeometry().Jacobian(J0, this_integration_method);
  MathUtils<double>::InvertMatrix(J0[PointNumber], invJ0[PointNumber],
                                  detJ0[PointNumber]);

  Matrix DN_DX = prod(DN_De[PointNumber], invJ0[PointNumber]);

  const SizeType num_nodes = rgeom.PointsNumber();

  MatrixType B;

  if (Dimension == 2) {
    B = ZeroMatrix(3, num_nodes * 2);

    SizeType index = 0;
    for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
      B(0, index + 0) = DN_DX(i_node, 0);
      B(0, index + 1) = 0.0;
      B(1, index + 0) = 0.0;
      B(1, index + 1) = DN_DX(i_node, 1);
      B(2, index + 0) = DN_DX(i_node, 1);
      B(2, index + 1) = DN_DX(i_node, 0);
      index += 2;
    }
  }

  else if (Dimension == 3) {
    B = ZeroMatrix(6, num_nodes * 3);

    SizeType index = 0;
    for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
      B(0, index + 0) = DN_DX(i_node, 0);
      B(1, index + 1) = DN_DX(i_node, 1);
      B(2, index + 2) = DN_DX(i_node, 2);
      B(3, index + 0) = DN_DX(i_node, 1);
      B(3, index + 1) = DN_DX(i_node, 0);
      B(4, index + 1) = DN_DX(i_node, 2);
      B(4, index + 2) = DN_DX(i_node, 1);
      B(5, index + 0) = DN_DX(i_node, 2);
      B(5, index + 2) = DN_DX(i_node, 0);
      index += 3;
    }
  }

  return B;

  KRATOS_CATCH("");
}

} // Namespace Kratos
