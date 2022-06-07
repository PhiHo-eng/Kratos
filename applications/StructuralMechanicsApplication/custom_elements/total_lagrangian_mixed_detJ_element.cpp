// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Riccardo Rossi
//                   Guglielmo Scovazzi
//

// System includes

// External includes

// Project includes
#include "utilities/element_size_calculator.h"
#include "utilities/geometry_utilities.h"
#include "utilities/math_utils.h"
#include "includes/checks.h"

// Application includes
#include "custom_elements/total_lagrangian_mixed_detJ_element.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"

namespace Kratos
{

template<std::size_t TDim>
Element::Pointer TotalLagrangianMixedDetJElement<TDim>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<TotalLagrangianMixedDetJElement<TDim>>(NewId, GetGeometry().Create(ThisNodes), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
Element::Pointer TotalLagrangianMixedDetJElement<TDim>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<TotalLagrangianMixedDetJElement<TDim>>(NewId, pGeom, pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
Element::Pointer TotalLagrangianMixedDetJElement<TDim>::Clone(
    IndexType NewId,
    NodesArrayType const& rThisNodes) const
{
    KRATOS_TRY

    TotalLagrangianMixedDetJElement::Pointer p_new_elem = Kratos::make_intrusive<TotalLagrangianMixedDetJElement>(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    // Currently selected integration methods
    p_new_elem->SetIntegrationMethod(mThisIntegrationMethod);

    // The vector containing the constitutive laws
    p_new_elem->SetConstitutiveLawVector(mConstitutiveLawVector);

    return p_new_elem;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void TotalLagrangianMixedDetJElement<2>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    if (rResult.size() != LocalSize){
        rResult.resize(LocalSize);
    }

    const auto& r_geometry = GetGeometry();
    const IndexType disp_pos = r_geometry[0].GetDofPosition(DISPLACEMENT_X);
    const IndexType det_F_pos = r_geometry[0].GetDofPosition(DETERMINANT_F);

    IndexType aux_index = 0;
    for (IndexType i_node = 0; i_node < NumNodes; ++i_node) {
        rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_X, disp_pos).EquationId();
        rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_Y, disp_pos + 1).EquationId();
        rResult[aux_index++] = r_geometry[i_node].GetDof(DETERMINANT_F, det_F_pos).EquationId();
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void TotalLagrangianMixedDetJElement<3>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    if (rResult.size() != LocalSize){
        rResult.resize(LocalSize);
    }

    const auto& r_geometry = GetGeometry();
    const IndexType disp_pos = r_geometry[0].GetDofPosition(DISPLACEMENT_X);
    const IndexType det_F_pos = r_geometry[0].GetDofPosition(DETERMINANT_F);

    IndexType aux_index = 0;
    for (IndexType i_node = 0; i_node < NumNodes; ++i_node) {
        rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_X, disp_pos).EquationId();
        rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_Y, disp_pos + 1).EquationId();
        rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_Z, disp_pos + 2).EquationId();
        rResult[aux_index++] = r_geometry[i_node].GetDof(DETERMINANT_F, det_F_pos).EquationId();
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void TotalLagrangianMixedDetJElement<2>::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    if (rElementalDofList.size() != LocalSize){
        rElementalDofList.resize(LocalSize);
    }

    const auto& r_geometry = GetGeometry();
    for(IndexType i = 0; i < NumNodes; ++i) {
        rElementalDofList[i * BlockSize] = r_geometry[i].pGetDof(DISPLACEMENT_X);
        rElementalDofList[i * BlockSize + 1] = r_geometry[i].pGetDof(DISPLACEMENT_Y);
        rElementalDofList[i * BlockSize + 2] = r_geometry[i].pGetDof(DETERMINANT_F);
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void TotalLagrangianMixedDetJElement<3>::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    if (rElementalDofList.size() != LocalSize){
        rElementalDofList.resize(LocalSize);
    }

    const auto& r_geometry = GetGeometry();
    for(IndexType i = 0; i < NumNodes; ++i){
        rElementalDofList[i * BlockSize] = r_geometry[i].pGetDof(DISPLACEMENT_X);
        rElementalDofList[i * BlockSize + 1] = r_geometry[i].pGetDof(DISPLACEMENT_Y);
        rElementalDofList[i * BlockSize + 2] = r_geometry[i].pGetDof(DISPLACEMENT_Z);
        rElementalDofList[i * BlockSize + 3] = r_geometry[i].pGetDof(DETERMINANT_F);
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void TotalLagrangianMixedDetJElement<TDim>::Initialize(const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY

    // Initialization should not be done again in a restart!
    if (!rCurrentProcessInfo[IS_RESTARTED]) {
        // Integration method initialization
        mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_2;
        const auto& r_integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

        // Constitutive Law Vector initialisation
        if (mConstitutiveLawVector.size() != r_integration_points.size()) {
            mConstitutiveLawVector.resize(r_integration_points.size());
        }

        // Initialize material
        InitializeMaterial();
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void TotalLagrangianMixedDetJElement<TDim>::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();
    const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());

    // Create the kinematics container and fill the nodal data
    KinematicVariables kinematic_variables;
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
        const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        for (IndexType d = 0; d < dim; ++d) {
            kinematic_variables.Displacements(i_node, d) = r_disp[d];
        }
        kinematic_variables.JacobianDeterminant[i_node] = r_geometry[i_node].FastGetSolutionStepValue(DETERMINANT_F);
    }

    // Set te constitutive law values
    ConstitutiveVariables constitutive_variables;
    ConstitutiveLaw::Parameters cons_law_values(GetGeometry(), GetProperties(), rCurrentProcessInfo);
    auto& r_cons_law_options = cons_law_values.GetOptions();
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // Call the initialize material response
    for (IndexType i_gauss = 0; i_gauss < mConstitutiveLawVector.size(); ++i_gauss) {
        // Recompute the kinematics
        CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());

        // Set the constitutive variables
        SetConstitutiveVariables(kinematic_variables, constitutive_variables, cons_law_values, i_gauss, r_integration_points);

        // Call the constitutive law to update material variables
        mConstitutiveLawVector[i_gauss]->InitializeMaterialResponseCauchy(cons_law_values);
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void TotalLagrangianMixedDetJElement<TDim>::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();
    const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());

    // Create the kinematics container and fill the nodal data
    KinematicVariables kinematic_variables;
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
        const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        for (IndexType d = 0; d < dim; ++d) {
            kinematic_variables.Displacements(i_node, d) = r_disp[d];
        }
        kinematic_variables.JacobianDeterminant[i_node] = r_geometry[i_node].FastGetSolutionStepValue(DETERMINANT_F);
    }

    // Set te constitutive law values
    ConstitutiveVariables constitutive_variables;
    ConstitutiveLaw::Parameters cons_law_values(r_geometry, GetProperties(), rCurrentProcessInfo);
    auto& r_cons_law_options = cons_law_values.GetOptions();
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // Call the initialize material response
    for (IndexType i_gauss = 0; i_gauss < mConstitutiveLawVector.size(); ++i_gauss) {
        // Recompute the kinematics
        CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());

        // Set the constitutive variables
        SetConstitutiveVariables(kinematic_variables, constitutive_variables, cons_law_values, i_gauss, r_integration_points);

        // Call the constitutive law to update material variables
        mConstitutiveLawVector[i_gauss]->FinalizeMaterialResponseCauchy(cons_law_values);
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void TotalLagrangianMixedDetJElement<2>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geometry = GetGeometry();
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType block_size = dim + 1;
    const SizeType matrix_size = block_size * n_nodes;
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();

    // Check RHS size
    if (rRightHandSideVector.size() != matrix_size) {
        rRightHandSideVector.resize(matrix_size, false);
    }

    // Check LHS size
    if (rLeftHandSideMatrix.size1() != matrix_size || rLeftHandSideMatrix.size2() != matrix_size) {
        rLeftHandSideMatrix.resize(matrix_size, matrix_size, false);
    }

    // Create the kinematics container and fill the nodal data
    KinematicVariables kinematic_variables;
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
        const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        for (IndexType d = 0; d < dim; ++d) {
            kinematic_variables.Displacements(i_node * dim + d) = r_disp[d];
        }
        kinematic_variables.JacobianDeterminant[i_node] = r_geometry[i_node].FastGetSolutionStepValue(DETERMINANT_F);
    }

    // Create the constitutive variables and values containers
    ConstitutiveVariables constitutive_variables;
    ConstitutiveLaw::Parameters cons_law_values(r_geometry, GetProperties(), rCurrentProcessInfo);
    auto& r_cons_law_options = cons_law_values.GetOptions();
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // Calculate the RHS and LHS contributions
    rLeftHandSideMatrix.clear();
    rRightHandSideVector.clear();

    // Calculate stabilization constant
    const double c_tau = 2.0;
    const double h = ElementSizeCalculator<2,NumNodes>::AverageElementSize(r_geometry);
    const double mu = 1.0; //FIXME: This is the Lame constant. Compute it.
    const double tau = c_tau * std::pow(h,2) / (2.0 * mu);

    // Set data for the body force calculation
    BoundedMatrix<double, NumNodes, 2> b;
    for (IndexType i_node = 0; i_node < NumNodes; ++i_node) {
        const array_1d<double,3>& r_b_i = r_geometry[i_node].FastGetSolutionStepValue(BODY_FORCE);
        for (IndexType d = 0; d < 2; ++d) {
            b(i_node, d) = r_b_i[d];
        }
    }
    const double rho0 = GetProperties().GetValue(DENSITY);

    // Set the auxiliary references matching the automatic differentiation symbols
    const auto& N = kinematic_variables.N;
    const auto& DN = kinematic_variables.DN_DX;
    const auto& u = kinematic_variables.Displacements;
    const auto& th = kinematic_variables.JacobianDeterminant;
    const auto& S = constitutive_variables.StressVector;
    const auto& C = constitutive_variables.ConstitutiveMatrix;

    // Aux RHS and LHS
    //TODO: To be removed
    BoundedVector<double,LocalSize> rhs;
    BoundedMatrix<double,LocalSize,LocalSize> lhs;

    const SizeType n_gauss = r_geometry.IntegrationPointsNumber(GetIntegrationMethod());
    const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());
    for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
        // Calculate kinematics
        CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());
        const double w_gauss = kinematic_variables.detJ0 * r_integration_points[i_gauss].Weight();

        // Calculate the constitutive response
        CalculateConstitutiveVariables(kinematic_variables, constitutive_variables, cons_law_values, i_gauss, r_geometry.IntegrationPoints(this->GetIntegrationMethod()), ConstitutiveLaw::StressMeasure_Cauchy);

        // Calculate and add the LHS Gauss point contributions
        const double clhs0 = DN(0,1)*u(0,0);
const double clhs1 = DN(1,1)*u(1,0);
const double clhs2 = DN(2,1)*u(2,0);
const double clhs3 = clhs0 + clhs1 + clhs2;
const double clhs4 = DN(0,1)*clhs3;
const double clhs5 = C(0,1)*clhs4;
const double clhs6 = DN(0,0)*u(0,0);
const double clhs7 = DN(1,0)*u(1,0);
const double clhs8 = DN(2,0)*u(2,0);
const double clhs9 = clhs6 + clhs7 + clhs8 + 1;
const double clhs10 = DN(0,0)*clhs9;
const double clhs11 = C(0,0)*clhs10;
const double clhs12 = DN(0,0)*clhs3 + DN(0,1)*clhs9;
const double clhs13 = C(0,2)*clhs12;
const double clhs14 = 1.0*clhs11 + clhs13 + 1.0*clhs5;
const double clhs15 = N[0]*th[0] + N[1]*th[1] + N[2]*th[2];
const double clhs16 = 1.0/clhs15;
const double clhs17 = clhs16*clhs9;
const double clhs18 = C(1,2)*clhs4;
const double clhs19 = C(0,2)*clhs10;
const double clhs20 = C(2,2)*clhs12;
const double clhs21 = 1.0*clhs18 + 1.0*clhs19 + clhs20;
const double clhs22 = clhs16*clhs3;
const double clhs23 = DN(0,0)*u(0,1);
const double clhs24 = DN(1,0)*u(1,1);
const double clhs25 = DN(2,0)*u(2,1);
const double clhs26 = clhs23 + clhs24 + clhs25;
const double clhs27 = pow(clhs26, 2);
const double clhs28 = pow(clhs9, 2);
const double clhs29 = 0.5*clhs16*clhs27 + 0.5*clhs16*clhs28 - 0.5;
const double clhs30 = pow(clhs3, 2);
const double clhs31 = DN(0,1)*u(0,1);
const double clhs32 = DN(1,1)*u(1,1);
const double clhs33 = DN(2,1)*u(2,1);
const double clhs34 = clhs31 + clhs32 + clhs33;
const double clhs35 = clhs34 + 1;
const double clhs36 = pow(clhs35, 2);
const double clhs37 = 0.5*clhs16*clhs30 + 0.5*clhs16*clhs36 - 0.5;
const double clhs38 = clhs26*clhs35 + clhs3*clhs9;
const double clhs39 = clhs16*clhs38;
const double clhs40 = C(0,0)*clhs29 + C(0,1)*clhs37 + C(0,2)*clhs39;
const double clhs41 = C(0,2)*clhs29 + C(1,2)*clhs37 + C(2,2)*clhs39;
const double clhs42 = DN(0,0)*clhs40 + DN(0,1)*clhs41;
const double clhs43 = clhs14*clhs17 + clhs21*clhs22 + clhs42;
const double clhs44 = C(1,1)*clhs4;
const double clhs45 = C(0,1)*clhs10;
const double clhs46 = C(1,2)*clhs12;
const double clhs47 = 1.0*clhs44 + 1.0*clhs45 + clhs46;
const double clhs48 = C(0,1)*clhs29 + C(1,1)*clhs37 + C(1,2)*clhs39;
const double clhs49 = DN(0,0)*clhs41 + DN(0,1)*clhs48;
const double clhs50 = clhs17*clhs21 + clhs22*clhs47 + clhs49;
const double clhs51 = DN(0,0)*clhs26;
const double clhs52 = C(0,2)*clhs51;
const double clhs53 = DN(0,1)*clhs35;
const double clhs54 = C(1,2)*clhs53;
const double clhs55 = DN(0,0)*clhs35 + DN(0,1)*clhs26;
const double clhs56 = C(2,2)*clhs55;
const double clhs57 = 1.0*clhs52 + 1.0*clhs54 + clhs56;
const double clhs58 = C(0,0)*clhs51;
const double clhs59 = C(0,1)*clhs53;
const double clhs60 = C(0,2)*clhs55;
const double clhs61 = 1.0*clhs58 + 1.0*clhs59 + clhs60;
const double clhs62 = clhs3*clhs57 + clhs61*clhs9;
const double clhs63 = C(0,1)*clhs51;
const double clhs64 = C(1,1)*clhs53;
const double clhs65 = C(1,2)*clhs55;
const double clhs66 = 1.0*clhs63 + 1.0*clhs64 + clhs65;
const double clhs67 = clhs3*clhs66 + clhs57*clhs9;
const double clhs68 = pow(clhs15, -2.0);
const double clhs69 = clhs27 + clhs28;
const double clhs70 = C(0,2)*clhs69;
const double clhs71 = clhs30 + clhs36;
const double clhs72 = C(1,2)*clhs71;
const double clhs73 = 1.0*clhs38;
const double clhs74 = C(2,2)*clhs73 + 0.5*clhs70 + 0.5*clhs72;
const double clhs75 = C(0,0)*clhs69;
const double clhs76 = C(0,1)*clhs71;
const double clhs77 = C(0,2)*clhs73 + 0.5*clhs75 + 0.5*clhs76;
const double clhs78 = clhs3*clhs74 + clhs77*clhs9;
const double clhs79 = C(0,1)*clhs69;
const double clhs80 = C(1,1)*clhs71;
const double clhs81 = C(1,2)*clhs73 + 0.5*clhs79 + 0.5*clhs80;
const double clhs82 = clhs3*clhs81 + clhs74*clhs9;
const double clhs83 = clhs68*(DN(0,0)*clhs78 + DN(0,1)*clhs82);
const double clhs84 = DN(1,1)*clhs3;
const double clhs85 = C(0,1)*clhs84;
const double clhs86 = DN(1,0)*clhs9;
const double clhs87 = C(0,0)*clhs86;
const double clhs88 = DN(1,0)*clhs3 + DN(1,1)*clhs9;
const double clhs89 = C(0,2)*clhs88;
const double clhs90 = 1.0*clhs85 + 1.0*clhs87 + clhs89;
const double clhs91 = C(1,2)*clhs84;
const double clhs92 = C(0,2)*clhs86;
const double clhs93 = C(2,2)*clhs88;
const double clhs94 = 1.0*clhs91 + 1.0*clhs92 + clhs93;
const double clhs95 = DN(1,0)*clhs40 + DN(1,1)*clhs41;
const double clhs96 = clhs17*clhs90 + clhs22*clhs94 + clhs95;
const double clhs97 = C(1,1)*clhs84;
const double clhs98 = C(0,1)*clhs86;
const double clhs99 = C(1,2)*clhs88;
const double clhs100 = 1.0*clhs97 + 1.0*clhs98 + clhs99;
const double clhs101 = DN(1,0)*clhs41 + DN(1,1)*clhs48;
const double clhs102 = clhs100*clhs22 + clhs101 + clhs17*clhs94;
const double clhs103 = DN(1,0)*clhs26;
const double clhs104 = C(0,2)*clhs103;
const double clhs105 = DN(1,1)*clhs35;
const double clhs106 = C(1,2)*clhs105;
const double clhs107 = DN(1,0)*clhs35 + DN(1,1)*clhs26;
const double clhs108 = C(2,2)*clhs107;
const double clhs109 = 1.0*clhs104 + 1.0*clhs106 + clhs108;
const double clhs110 = C(0,0)*clhs103;
const double clhs111 = C(0,1)*clhs105;
const double clhs112 = C(0,2)*clhs107;
const double clhs113 = 1.0*clhs110 + 1.0*clhs111 + clhs112;
const double clhs114 = clhs109*clhs3 + clhs113*clhs9;
const double clhs115 = C(0,1)*clhs103;
const double clhs116 = C(1,1)*clhs105;
const double clhs117 = C(1,2)*clhs107;
const double clhs118 = 1.0*clhs115 + 1.0*clhs116 + clhs117;
const double clhs119 = clhs109*clhs9 + clhs118*clhs3;
const double clhs120 = DN(2,1)*clhs3;
const double clhs121 = C(0,1)*clhs120;
const double clhs122 = DN(2,0)*clhs9;
const double clhs123 = C(0,0)*clhs122;
const double clhs124 = DN(2,0)*clhs3 + DN(2,1)*clhs9;
const double clhs125 = C(0,2)*clhs124;
const double clhs126 = 1.0*clhs121 + 1.0*clhs123 + clhs125;
const double clhs127 = C(1,2)*clhs120;
const double clhs128 = C(0,2)*clhs122;
const double clhs129 = C(2,2)*clhs124;
const double clhs130 = 1.0*clhs127 + 1.0*clhs128 + clhs129;
const double clhs131 = DN(2,0)*clhs40 + DN(2,1)*clhs41;
const double clhs132 = clhs126*clhs17 + clhs130*clhs22 + clhs131;
const double clhs133 = C(1,1)*clhs120;
const double clhs134 = C(0,1)*clhs122;
const double clhs135 = C(1,2)*clhs124;
const double clhs136 = 1.0*clhs133 + 1.0*clhs134 + clhs135;
const double clhs137 = DN(2,0)*clhs41 + DN(2,1)*clhs48;
const double clhs138 = clhs130*clhs17 + clhs136*clhs22 + clhs137;
const double clhs139 = DN(2,0)*clhs26;
const double clhs140 = C(0,2)*clhs139;
const double clhs141 = DN(2,1)*clhs35;
const double clhs142 = C(1,2)*clhs141;
const double clhs143 = DN(2,0)*clhs35 + DN(2,1)*clhs26;
const double clhs144 = C(2,2)*clhs143;
const double clhs145 = 1.0*clhs140 + 1.0*clhs142 + clhs144;
const double clhs146 = C(0,0)*clhs139;
const double clhs147 = C(0,1)*clhs141;
const double clhs148 = C(0,2)*clhs143;
const double clhs149 = 1.0*clhs146 + 1.0*clhs147 + clhs148;
const double clhs150 = clhs145*clhs3 + clhs149*clhs9;
const double clhs151 = C(0,1)*clhs139;
const double clhs152 = C(1,1)*clhs141;
const double clhs153 = C(1,2)*clhs143;
const double clhs154 = 1.0*clhs151 + 1.0*clhs152 + clhs153;
const double clhs155 = clhs145*clhs9 + clhs154*clhs3;
const double clhs156 = clhs14*clhs26 + clhs21*clhs35;
const double clhs157 = clhs21*clhs26 + clhs35*clhs47;
const double clhs158 = clhs16*clhs26;
const double clhs159 = clhs16*clhs35;
const double clhs160 = clhs158*clhs61 + clhs159*clhs57 + clhs42;
const double clhs161 = clhs158*clhs57 + clhs159*clhs66 + clhs49;
const double clhs162 = clhs26*clhs77 + clhs35*clhs74;
const double clhs163 = clhs26*clhs74 + clhs35*clhs81;
const double clhs164 = clhs68*(DN(0,0)*clhs162 + DN(0,1)*clhs163);
const double clhs165 = clhs26*clhs90 + clhs35*clhs94;
const double clhs166 = clhs100*clhs35 + clhs26*clhs94;
const double clhs167 = clhs109*clhs159 + clhs113*clhs158 + clhs95;
const double clhs168 = clhs101 + clhs109*clhs158 + clhs118*clhs159;
const double clhs169 = clhs126*clhs26 + clhs130*clhs35;
const double clhs170 = clhs130*clhs26 + clhs136*clhs35;
const double clhs171 = clhs131 + clhs145*clhs159 + clhs149*clhs158;
const double clhs172 = clhs137 + clhs145*clhs158 + clhs154*clhs159;
const double clhs173 = DN(0,0)*clhs32 + DN(0,0)*clhs33 + DN(0,0) - DN(0,1)*clhs24 - DN(0,1)*clhs25;
const double clhs174 = clhs11 + clhs13 + clhs5;
const double clhs175 = clhs18 + clhs19 + clhs20;
const double clhs176 = 2*clhs38;
const double clhs177 = C(0,2)*clhs176 + clhs75 + clhs76;
const double clhs178 = DN(0,0)*clhs177;
const double clhs179 = -clhs0*clhs24 - clhs0*clhs25 - clhs1*clhs23 - clhs1*clhs25 - clhs2*clhs23 - clhs2*clhs24 + clhs31*clhs7 + clhs31*clhs8 + clhs32*clhs6 + clhs32*clhs8 + clhs33*clhs6 + clhs33*clhs7 + clhs34 + clhs9;
const double clhs180 = 0.25/clhs179;
const double clhs181 = clhs173*clhs180;
const double clhs182 = C(2,2)*clhs176 + clhs70 + clhs72;
const double clhs183 = DN(0,1)*clhs182;
const double clhs184 = DN(0,0)*th[0] + DN(1,0)*th[1] + DN(2,0)*th[2];
const double clhs185 = sqrt(clhs179/clhs15);
const double clhs186 = clhs185*tau;
const double clhs187 = clhs16*clhs186;
const double clhs188 = clhs184*clhs187;
const double clhs189 = clhs44 + clhs45 + clhs46;
const double clhs190 = DN(0,0)*clhs182;
const double clhs191 = C(1,2)*clhs176 + clhs79 + clhs80;
const double clhs192 = DN(0,1)*clhs191;
const double clhs193 = DN(0,1)*th[0] + DN(1,1)*th[1] + DN(2,1)*th[2];
const double clhs194 = clhs187*clhs193;
const double clhs195 = -DN(0,0)*clhs1 - DN(0,0)*clhs2 + DN(0,1)*clhs7 + DN(0,1)*clhs8 + DN(0,1);
const double clhs196 = clhs58 + clhs59 + clhs60;
const double clhs197 = clhs52 + clhs54 + clhs56;
const double clhs198 = clhs180*clhs195;
const double clhs199 = clhs63 + clhs64 + clhs65;
const double clhs200 = (1.0/2.0)*clhs187;
const double clhs201 = clhs200*(clhs178 + clhs183);
const double clhs202 = clhs200*(clhs190 + clhs192);
const double clhs203 = (1.0/2.0)*clhs186*clhs68;
const double clhs204 = N[0]*clhs203;
const double clhs205 = 2.0*clhs38;
const double clhs206 = C(0,2)*clhs205 + 1.0*clhs75 + 1.0*clhs76;
const double clhs207 = C(2,2)*clhs205 + 1.0*clhs70 + 1.0*clhs72;
const double clhs208 = clhs184*(DN(0,0)*clhs206 + DN(0,1)*clhs207 + 0.5*clhs178 + 0.5*clhs183);
const double clhs209 = C(1,2)*clhs205 + 1.0*clhs79 + 1.0*clhs80;
const double clhs210 = clhs193*(DN(0,0)*clhs207 + DN(0,1)*clhs209 + 0.5*clhs190 + 0.5*clhs192);
const double clhs211 = N[0]*b(0,1) + N[1]*b(1,1) + N[2]*b(2,1);
const double clhs212 = rho0*tau;
const double clhs213 = clhs212*(DN(0,0)*DN(1,1) - DN(0,1)*DN(1,0));
const double clhs214 = clhs211*clhs213;
const double clhs215 = DN(1,0)*clhs31 + DN(1,0)*clhs33 + DN(1,0) - DN(1,1)*clhs23 - DN(1,1)*clhs25;
const double clhs216 = clhs85 + clhs87 + clhs89;
const double clhs217 = clhs91 + clhs92 + clhs93;
const double clhs218 = clhs180*clhs215;
const double clhs219 = clhs97 + clhs98 + clhs99;
const double clhs220 = N[0]*b(0,0) + N[1]*b(1,0) + N[2]*b(2,0);
const double clhs221 = clhs213*clhs220;
const double clhs222 = -DN(1,0)*clhs0 - DN(1,0)*clhs2 + DN(1,1)*clhs6 + DN(1,1)*clhs8 + DN(1,1);
const double clhs223 = clhs110 + clhs111 + clhs112;
const double clhs224 = clhs104 + clhs106 + clhs108;
const double clhs225 = clhs180*clhs222;
const double clhs226 = clhs115 + clhs116 + clhs117;
const double clhs227 = N[0]*N[1];
const double clhs228 = N[1]*clhs203;
const double clhs229 = clhs212*(DN(0,0)*DN(2,1) - DN(0,1)*DN(2,0));
const double clhs230 = clhs211*clhs229;
const double clhs231 = DN(2,0)*clhs31 + DN(2,0)*clhs32 + DN(2,0) - DN(2,1)*clhs23 - DN(2,1)*clhs24;
const double clhs232 = clhs121 + clhs123 + clhs125;
const double clhs233 = clhs127 + clhs128 + clhs129;
const double clhs234 = clhs180*clhs231;
const double clhs235 = clhs133 + clhs134 + clhs135;
const double clhs236 = clhs220*clhs229;
const double clhs237 = -DN(2,0)*clhs0 - DN(2,0)*clhs1 + DN(2,1)*clhs6 + DN(2,1)*clhs7 + DN(2,1);
const double clhs238 = clhs146 + clhs147 + clhs148;
const double clhs239 = clhs140 + clhs142 + clhs144;
const double clhs240 = clhs180*clhs237;
const double clhs241 = clhs151 + clhs152 + clhs153;
const double clhs242 = N[0]*N[2];
const double clhs243 = N[2]*clhs203;
const double clhs244 = clhs68*(DN(1,0)*clhs78 + DN(1,1)*clhs82);
const double clhs245 = clhs68*(DN(1,0)*clhs162 + DN(1,1)*clhs163);
const double clhs246 = DN(1,0)*clhs177;
const double clhs247 = DN(1,1)*clhs182;
const double clhs248 = DN(1,0)*clhs182;
const double clhs249 = DN(1,1)*clhs191;
const double clhs250 = clhs200*(clhs246 + clhs247);
const double clhs251 = clhs200*(clhs248 + clhs249);
const double clhs252 = clhs184*(DN(1,0)*clhs206 + DN(1,1)*clhs207 + 0.5*clhs246 + 0.5*clhs247);
const double clhs253 = clhs193*(DN(1,0)*clhs207 + DN(1,1)*clhs209 + 0.5*clhs248 + 0.5*clhs249);
const double clhs254 = clhs212*(DN(1,0)*DN(2,1) - DN(1,1)*DN(2,0));
const double clhs255 = clhs211*clhs254;
const double clhs256 = clhs220*clhs254;
const double clhs257 = N[1]*N[2];
const double clhs258 = clhs68*(DN(2,0)*clhs78 + DN(2,1)*clhs82);
const double clhs259 = clhs68*(DN(2,0)*clhs162 + DN(2,1)*clhs163);
const double clhs260 = DN(2,0)*clhs177;
const double clhs261 = DN(2,1)*clhs182;
const double clhs262 = DN(2,0)*clhs182;
const double clhs263 = DN(2,1)*clhs191;
const double clhs264 = clhs200*(clhs260 + clhs261);
const double clhs265 = clhs200*(clhs262 + clhs263);
const double clhs266 = clhs184*(DN(2,0)*clhs206 + DN(2,1)*clhs207 + 0.5*clhs260 + 0.5*clhs261);
const double clhs267 = clhs193*(DN(2,0)*clhs207 + DN(2,1)*clhs209 + 0.5*clhs262 + 0.5*clhs263);
lhs(0,0)=-DN(0,0)*clhs43 - DN(0,1)*clhs50;
lhs(0,1)=-clhs16*(DN(0,0)*clhs62 + DN(0,1)*clhs67);
lhs(0,2)=N[0]*clhs83;
lhs(0,3)=-DN(0,0)*clhs96 - DN(0,1)*clhs102;
lhs(0,4)=-clhs16*(DN(0,0)*clhs114 + DN(0,1)*clhs119);
lhs(0,5)=N[1]*clhs83;
lhs(0,6)=-DN(0,0)*clhs132 - DN(0,1)*clhs138;
lhs(0,7)=-clhs16*(DN(0,0)*clhs150 + DN(0,1)*clhs155);
lhs(0,8)=N[2]*clhs83;
lhs(1,0)=-clhs16*(DN(0,0)*clhs156 + DN(0,1)*clhs157);
lhs(1,1)=-DN(0,0)*clhs160 - DN(0,1)*clhs161;
lhs(1,2)=N[0]*clhs164;
lhs(1,3)=-clhs16*(DN(0,0)*clhs165 + DN(0,1)*clhs166);
lhs(1,4)=-DN(0,0)*clhs167 - DN(0,1)*clhs168;
lhs(1,5)=N[1]*clhs164;
lhs(1,6)=-clhs16*(DN(0,0)*clhs169 + DN(0,1)*clhs170);
lhs(1,7)=-DN(0,0)*clhs171 - DN(0,1)*clhs172;
lhs(1,8)=N[2]*clhs164;
lhs(2,0)=-N[0]*clhs173 + clhs188*(DN(0,0)*clhs174 + DN(0,1)*clhs175 + clhs178*clhs181 + clhs181*clhs183) + clhs194*(DN(0,0)*clhs175 + DN(0,1)*clhs189 + clhs181*clhs190 + clhs181*clhs192);
lhs(2,1)=-N[0]*clhs195 + clhs188*(DN(0,0)*clhs196 + DN(0,1)*clhs197 + clhs178*clhs198 + clhs183*clhs198) + clhs194*(DN(0,0)*clhs197 + DN(0,1)*clhs199 + clhs190*clhs198 + clhs192*clhs198);
lhs(2,2)=DN(0,0)*clhs201 + DN(0,1)*clhs202 + pow(N[0], 2) - clhs204*clhs208 - clhs204*clhs210;
lhs(2,3)=-N[0]*clhs215 + clhs16*clhs184*clhs185*tau*(DN(0,0)*clhs216 + DN(0,1)*clhs217 + clhs178*clhs218 + clhs183*clhs218) + clhs16*clhs185*clhs193*tau*(DN(0,0)*clhs217 + DN(0,1)*clhs219 + clhs190*clhs218 + clhs192*clhs218) - clhs214;
lhs(2,4)=-N[0]*clhs222 + clhs188*(DN(0,0)*clhs223 + DN(0,1)*clhs224 + clhs178*clhs225 + clhs183*clhs225) + clhs194*(DN(0,0)*clhs224 + DN(0,1)*clhs226 + clhs190*clhs225 + clhs192*clhs225) + clhs221;
lhs(2,5)=DN(1,0)*clhs201 + DN(1,1)*clhs202 - clhs208*clhs228 - clhs210*clhs228 + clhs227;
lhs(2,6)=-N[0]*clhs231 + clhs16*clhs184*clhs185*tau*(DN(0,0)*clhs232 + DN(0,1)*clhs233 + clhs178*clhs234 + clhs183*clhs234) + clhs16*clhs185*clhs193*tau*(DN(0,0)*clhs233 + DN(0,1)*clhs235 + clhs190*clhs234 + clhs192*clhs234) - clhs230;
lhs(2,7)=-N[0]*clhs237 + clhs188*(DN(0,0)*clhs238 + DN(0,1)*clhs239 + clhs178*clhs240 + clhs183*clhs240) + clhs194*(DN(0,0)*clhs239 + DN(0,1)*clhs241 + clhs190*clhs240 + clhs192*clhs240) + clhs236;
lhs(2,8)=DN(2,0)*clhs201 + DN(2,1)*clhs202 - clhs208*clhs243 - clhs210*clhs243 + clhs242;
lhs(3,0)=-DN(1,0)*clhs43 - DN(1,1)*clhs50;
lhs(3,1)=-clhs16*(DN(1,0)*clhs62 + DN(1,1)*clhs67);
lhs(3,2)=N[0]*clhs244;
lhs(3,3)=-DN(1,0)*clhs96 - DN(1,1)*clhs102;
lhs(3,4)=-clhs16*(DN(1,0)*clhs114 + DN(1,1)*clhs119);
lhs(3,5)=N[1]*clhs244;
lhs(3,6)=-DN(1,0)*clhs132 - DN(1,1)*clhs138;
lhs(3,7)=-clhs16*(DN(1,0)*clhs150 + DN(1,1)*clhs155);
lhs(3,8)=N[2]*clhs244;
lhs(4,0)=-clhs16*(DN(1,0)*clhs156 + DN(1,1)*clhs157);
lhs(4,1)=-DN(1,0)*clhs160 - DN(1,1)*clhs161;
lhs(4,2)=N[0]*clhs245;
lhs(4,3)=-clhs16*(DN(1,0)*clhs165 + DN(1,1)*clhs166);
lhs(4,4)=-DN(1,0)*clhs167 - DN(1,1)*clhs168;
lhs(4,5)=N[1]*clhs245;
lhs(4,6)=-clhs16*(DN(1,0)*clhs169 + DN(1,1)*clhs170);
lhs(4,7)=-DN(1,0)*clhs171 - DN(1,1)*clhs172;
lhs(4,8)=N[2]*clhs245;
lhs(5,0)=-N[1]*clhs173 + clhs188*(DN(1,0)*clhs174 + DN(1,1)*clhs175 + clhs181*clhs246 + clhs181*clhs247) + clhs194*(DN(1,0)*clhs175 + DN(1,1)*clhs189 + clhs181*clhs248 + clhs181*clhs249) + clhs214;
lhs(5,1)=-N[1]*clhs195 + clhs16*clhs184*clhs185*tau*(DN(1,0)*clhs196 + DN(1,1)*clhs197 + clhs198*clhs246 + clhs198*clhs247) + clhs16*clhs185*clhs193*tau*(DN(1,0)*clhs197 + DN(1,1)*clhs199 + clhs198*clhs248 + clhs198*clhs249) - clhs221;
lhs(5,2)=DN(0,0)*clhs250 + DN(0,1)*clhs251 - clhs204*clhs252 - clhs204*clhs253 + clhs227;
lhs(5,3)=-N[1]*clhs215 + clhs188*(DN(1,0)*clhs216 + DN(1,1)*clhs217 + clhs218*clhs246 + clhs218*clhs247) + clhs194*(DN(1,0)*clhs217 + DN(1,1)*clhs219 + clhs218*clhs248 + clhs218*clhs249);
lhs(5,4)=-N[1]*clhs222 + clhs188*(DN(1,0)*clhs223 + DN(1,1)*clhs224 + clhs225*clhs246 + clhs225*clhs247) + clhs194*(DN(1,0)*clhs224 + DN(1,1)*clhs226 + clhs225*clhs248 + clhs225*clhs249);
lhs(5,5)=DN(1,0)*clhs250 + DN(1,1)*clhs251 + pow(N[1], 2) - clhs228*clhs252 - clhs228*clhs253;
lhs(5,6)=-N[1]*clhs231 + clhs16*clhs184*clhs185*tau*(DN(1,0)*clhs232 + DN(1,1)*clhs233 + clhs234*clhs246 + clhs234*clhs247) + clhs16*clhs185*clhs193*tau*(DN(1,0)*clhs233 + DN(1,1)*clhs235 + clhs234*clhs248 + clhs234*clhs249) - clhs255;
lhs(5,7)=-N[1]*clhs237 + clhs188*(DN(1,0)*clhs238 + DN(1,1)*clhs239 + clhs240*clhs246 + clhs240*clhs247) + clhs194*(DN(1,0)*clhs239 + DN(1,1)*clhs241 + clhs240*clhs248 + clhs240*clhs249) + clhs256;
lhs(5,8)=DN(2,0)*clhs250 + DN(2,1)*clhs251 - clhs243*clhs252 - clhs243*clhs253 + clhs257;
lhs(6,0)=-DN(2,0)*clhs43 - DN(2,1)*clhs50;
lhs(6,1)=-clhs16*(DN(2,0)*clhs62 + DN(2,1)*clhs67);
lhs(6,2)=N[0]*clhs258;
lhs(6,3)=-DN(2,0)*clhs96 - DN(2,1)*clhs102;
lhs(6,4)=-clhs16*(DN(2,0)*clhs114 + DN(2,1)*clhs119);
lhs(6,5)=N[1]*clhs258;
lhs(6,6)=-DN(2,0)*clhs132 - DN(2,1)*clhs138;
lhs(6,7)=-clhs16*(DN(2,0)*clhs150 + DN(2,1)*clhs155);
lhs(6,8)=N[2]*clhs258;
lhs(7,0)=-clhs16*(DN(2,0)*clhs156 + DN(2,1)*clhs157);
lhs(7,1)=-DN(2,0)*clhs160 - DN(2,1)*clhs161;
lhs(7,2)=N[0]*clhs259;
lhs(7,3)=-clhs16*(DN(2,0)*clhs165 + DN(2,1)*clhs166);
lhs(7,4)=-DN(2,0)*clhs167 - DN(2,1)*clhs168;
lhs(7,5)=N[1]*clhs259;
lhs(7,6)=-clhs16*(DN(2,0)*clhs169 + DN(2,1)*clhs170);
lhs(7,7)=-DN(2,0)*clhs171 - DN(2,1)*clhs172;
lhs(7,8)=N[2]*clhs259;
lhs(8,0)=-N[2]*clhs173 + clhs188*(DN(2,0)*clhs174 + DN(2,1)*clhs175 + clhs181*clhs260 + clhs181*clhs261) + clhs194*(DN(2,0)*clhs175 + DN(2,1)*clhs189 + clhs181*clhs262 + clhs181*clhs263) + clhs230;
lhs(8,1)=-N[2]*clhs195 + clhs16*clhs184*clhs185*tau*(DN(2,0)*clhs196 + DN(2,1)*clhs197 + clhs198*clhs260 + clhs198*clhs261) + clhs16*clhs185*clhs193*tau*(DN(2,0)*clhs197 + DN(2,1)*clhs199 + clhs198*clhs262 + clhs198*clhs263) - clhs236;
lhs(8,2)=DN(0,0)*clhs264 + DN(0,1)*clhs265 - clhs204*clhs266 - clhs204*clhs267 + clhs242;
lhs(8,3)=-N[2]*clhs215 + clhs188*(DN(2,0)*clhs216 + DN(2,1)*clhs217 + clhs218*clhs260 + clhs218*clhs261) + clhs194*(DN(2,0)*clhs217 + DN(2,1)*clhs219 + clhs218*clhs262 + clhs218*clhs263) + clhs255;
lhs(8,4)=-N[2]*clhs222 + clhs16*clhs184*clhs185*tau*(DN(2,0)*clhs223 + DN(2,1)*clhs224 + clhs225*clhs260 + clhs225*clhs261) + clhs16*clhs185*clhs193*tau*(DN(2,0)*clhs224 + DN(2,1)*clhs226 + clhs225*clhs262 + clhs225*clhs263) - clhs256;
lhs(8,5)=DN(1,0)*clhs264 + DN(1,1)*clhs265 - clhs228*clhs266 - clhs228*clhs267 + clhs257;
lhs(8,6)=-N[2]*clhs231 + clhs188*(DN(2,0)*clhs232 + DN(2,1)*clhs233 + clhs234*clhs260 + clhs234*clhs261) + clhs194*(DN(2,0)*clhs233 + DN(2,1)*clhs235 + clhs234*clhs262 + clhs234*clhs263);
lhs(8,7)=-N[2]*clhs237 + clhs188*(DN(2,0)*clhs238 + DN(2,1)*clhs239 + clhs240*clhs260 + clhs240*clhs261) + clhs194*(DN(2,0)*clhs239 + DN(2,1)*clhs241 + clhs240*clhs262 + clhs240*clhs263);
lhs(8,8)=DN(2,0)*clhs264 + DN(2,1)*clhs265 + pow(N[2], 2) - clhs243*clhs266 - clhs243*clhs267;

        // Calculate and add the RHS Gauss point contribution
        const double crhs0 = N[0]*b(0,0) + N[1]*b(1,0) + N[2]*b(2,0);
const double crhs1 = N[0]*rho0;
const double crhs2 = DN(0,1)*u(0,0);
const double crhs3 = DN(1,1)*u(1,0);
const double crhs4 = DN(2,1)*u(2,0);
const double crhs5 = crhs2 + crhs3 + crhs4;
const double crhs6 = DN(0,0)*u(0,0);
const double crhs7 = DN(1,0)*u(1,0);
const double crhs8 = DN(2,0)*u(2,0);
const double crhs9 = crhs6 + crhs7 + crhs8 + 1;
const double crhs10 = S[0]*crhs9 + S[2]*crhs5;
const double crhs11 = S[1]*crhs5 + S[2]*crhs9;
const double crhs12 = N[0]*b(0,1) + N[1]*b(1,1) + N[2]*b(2,1);
const double crhs13 = DN(0,0)*u(0,1);
const double crhs14 = DN(1,0)*u(1,1);
const double crhs15 = DN(2,0)*u(2,1);
const double crhs16 = crhs13 + crhs14 + crhs15;
const double crhs17 = DN(0,1)*u(0,1);
const double crhs18 = DN(1,1)*u(1,1);
const double crhs19 = DN(2,1)*u(2,1);
const double crhs20 = crhs17 + crhs18 + crhs19;
const double crhs21 = crhs20 + 1;
const double crhs22 = S[0]*crhs16 + S[2]*crhs21;
const double crhs23 = S[1]*crhs21 + S[2]*crhs16;
const double crhs24 = rho0*tau;
const double crhs25 = crhs12*crhs24;
const double crhs26 = crhs0*crhs24;
const double crhs27 = N[0]*th[0];
const double crhs28 = N[1]*th[1];
const double crhs29 = N[2]*th[2];
const double crhs30 = -crhs13*crhs3 - crhs13*crhs4 - crhs14*crhs2 - crhs14*crhs4 - crhs15*crhs2 - crhs15*crhs3 + crhs17*crhs7 + crhs17*crhs8 + crhs18*crhs6 + crhs18*crhs8 + crhs19*crhs6 + crhs19*crhs7 + crhs20 + crhs9;
const double crhs31 = -crhs27 - crhs28 - crhs29 + crhs30;
const double crhs32 = pow(crhs16, 2) + pow(crhs9, 2);
const double crhs33 = pow(crhs21, 2) + pow(crhs5, 2);
const double crhs34 = 2*crhs16*crhs21 + 2*crhs5*crhs9;
const double crhs35 = C(0,0)*crhs32 + C(0,1)*crhs33 + C(0,2)*crhs34;
const double crhs36 = C(0,2)*crhs32 + C(1,2)*crhs33 + C(2,2)*crhs34;
const double crhs37 = crhs27 + crhs28 + crhs29;
const double crhs38 = (1.0/2.0)*1.0/crhs37*tau*sqrt(crhs30/crhs37);
const double crhs39 = crhs38*(DN(0,0)*th[0] + DN(1,0)*th[1] + DN(2,0)*th[2]);
const double crhs40 = C(0,1)*crhs32 + C(1,1)*crhs33 + C(1,2)*crhs34;
const double crhs41 = crhs38*(DN(0,1)*th[0] + DN(1,1)*th[1] + DN(2,1)*th[2]);
const double crhs42 = N[1]*rho0;
const double crhs43 = N[2]*rho0;
rhs[0]=DN(0,0)*crhs10 + DN(0,1)*crhs11 - crhs0*crhs1;
rhs[1]=DN(0,0)*crhs22 + DN(0,1)*crhs23 - crhs1*crhs12;
rhs[2]=N[0]*crhs31 - crhs25*(-DN(0,0)*crhs5 + DN(0,1)*crhs9) - crhs26*(DN(0,0)*crhs21 - DN(0,1)*crhs16) - crhs39*(DN(0,0)*crhs35 + DN(0,1)*crhs36) - crhs41*(DN(0,0)*crhs36 + DN(0,1)*crhs40);
rhs[3]=DN(1,0)*crhs10 + DN(1,1)*crhs11 - crhs0*crhs42;
rhs[4]=DN(1,0)*crhs22 + DN(1,1)*crhs23 - crhs12*crhs42;
rhs[5]=N[1]*crhs31 - crhs25*(-DN(1,0)*crhs5 + DN(1,1)*crhs9) - crhs26*(DN(1,0)*crhs21 - DN(1,1)*crhs16) - crhs39*(DN(1,0)*crhs35 + DN(1,1)*crhs36) - crhs41*(DN(1,0)*crhs36 + DN(1,1)*crhs40);
rhs[6]=DN(2,0)*crhs10 + DN(2,1)*crhs11 - crhs0*crhs43;
rhs[7]=DN(2,0)*crhs22 + DN(2,1)*crhs23 - crhs12*crhs43;
rhs[8]=N[2]*crhs31 - crhs25*(-DN(2,0)*crhs5 + DN(2,1)*crhs9) - crhs26*(DN(2,0)*crhs21 - DN(2,1)*crhs16) - crhs39*(DN(2,0)*crhs35 + DN(2,1)*crhs36) - crhs41*(DN(2,0)*crhs36 + DN(2,1)*crhs40);

        //TODO: Amend this once the assembly is done in the input arrays
        rLeftHandSideMatrix += w_gauss * lhs;
        rRightHandSideVector += w_gauss * rhs;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void TotalLagrangianMixedDetJElement<3>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geometry = GetGeometry();
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType block_size = dim + 1;
    const SizeType matrix_size = block_size * n_nodes;
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();

    // Check RHS size
    if (rRightHandSideVector.size() != matrix_size) {
        rRightHandSideVector.resize(matrix_size, false);
    }

    // Check LHS size
    if (rLeftHandSideMatrix.size1() != matrix_size || rLeftHandSideMatrix.size2() != matrix_size) {
        rLeftHandSideMatrix.resize(matrix_size, matrix_size, false);
    }

    // Create the kinematics container and fill the nodal data
    KinematicVariables kinematic_variables;
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
        const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        for (IndexType d = 0; d < dim; ++d) {
            kinematic_variables.Displacements(i_node * dim + d) = r_disp[d];
        }
        kinematic_variables.JacobianDeterminant[i_node] = r_geometry[i_node].FastGetSolutionStepValue(DETERMINANT_F);
    }

    // Create the constitutive variables and values containers
    ConstitutiveVariables constitutive_variables;
    ConstitutiveLaw::Parameters cons_law_values(r_geometry, GetProperties(), rCurrentProcessInfo);
    auto& r_cons_law_options = cons_law_values.GetOptions();
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // Calculate the RHS and LHS contributions
    rLeftHandSideMatrix.clear();
    rRightHandSideVector.clear();

    // Calculate stabilization constant
    const double c_tau = 2.0;
    const double h = ElementSizeCalculator<3,NumNodes>::AverageElementSize(r_geometry);
    const double mu = 1.0; //FIXME: This is the Lame constant. Compute it.
    const double tau = c_tau * std::pow(h,2) / (2.0 * mu);

    // Set data for the body force calculation
    BoundedMatrix<double, NumNodes, 3> b;
    for (IndexType i_node = 0; i_node < NumNodes; ++i_node) {
        const array_1d<double,3>& r_b_i = r_geometry[i_node].FastGetSolutionStepValue(BODY_FORCE);
        for (IndexType d = 0; d < 3; ++d) {
            b(i_node, d) = r_b_i[d];
        }
    }
    const double rho0 = GetProperties().GetValue(DENSITY);

    // Set the auxiliary references matching the automatic differentiation symbols
    const auto& N = kinematic_variables.N;
    const auto& DN = kinematic_variables.DN_DX;
    const auto& u = kinematic_variables.Displacements;
    const auto& th = kinematic_variables.JacobianDeterminant;
    const auto& S = constitutive_variables.StressVector;
    const auto& C = constitutive_variables.ConstitutiveMatrix;

    // Aux RHS and LHS
    //TODO: To be removed
    BoundedVector<double,LocalSize> rhs;
    BoundedMatrix<double,LocalSize,LocalSize> lhs;

    const SizeType n_gauss = r_geometry.IntegrationPointsNumber(GetIntegrationMethod());
    const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());
    for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
        // Calculate kinematics
        CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());
        const double w_gauss = kinematic_variables.detJ0 * r_integration_points[i_gauss].Weight();

        // Calculate the constitutive response
        CalculateConstitutiveVariables(kinematic_variables, constitutive_variables, cons_law_values, i_gauss, r_geometry.IntegrationPoints(this->GetIntegrationMethod()), ConstitutiveLaw::StressMeasure_Cauchy);

        // Calculate and add the LHS Gauss point contributions
        //substitute_lhs_3D_4N
        // Calculate and add the RHS Gauss point contribution
        //substitute_rhs_3D_4N
        //TODO: Amend this once the assembly is done in the input arrays
        rLeftHandSideMatrix += w_gauss * lhs;
        rRightHandSideVector += w_gauss * rhs;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void TotalLagrangianMixedDetJElement<2>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geometry = GetGeometry();
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType block_size = dim + 1;
    const SizeType matrix_size = block_size * n_nodes;
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();

    // Check LHS size
    if (rLeftHandSideMatrix.size1() != matrix_size || rLeftHandSideMatrix.size2() != matrix_size) {
        rLeftHandSideMatrix.resize(matrix_size, matrix_size, false);
    }

    // Create the kinematics container and fill the nodal data
    KinematicVariables kinematic_variables;
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
        const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        for (IndexType d = 0; d < dim; ++d) {
            kinematic_variables.Displacements(i_node * dim + d) = r_disp[d];
        }
        kinematic_variables.JacobianDeterminant[i_node] = r_geometry[i_node].FastGetSolutionStepValue(DETERMINANT_F);
    }

    // Create the constitutive variables and values containers
    ConstitutiveVariables constitutive_variables;
    ConstitutiveLaw::Parameters cons_law_values(r_geometry, GetProperties(), rCurrentProcessInfo);
    auto& r_cons_law_options = cons_law_values.GetOptions();
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // Calculate the RHS and LHS contributions
    rLeftHandSideMatrix.clear();

    // Calculate stabilization constant
    const double c_tau = 2.0;
    const double h = ElementSizeCalculator<2,NumNodes>::AverageElementSize(r_geometry);
    const double mu = 1.0; //FIXME: This is the Lame constant. Compute it.
    const double tau = c_tau * std::pow(h,2) / (2.0 * mu);

    // Set data for the body force calculation
    BoundedMatrix<double, NumNodes, 2> b;
    for (IndexType i_node = 0; i_node < NumNodes; ++i_node) {
        const array_1d<double,3>& r_b_i = r_geometry[i_node].FastGetSolutionStepValue(BODY_FORCE);
        for (IndexType d = 0; d < 2; ++d) {
            b(i_node, d) = r_b_i[d];
        }
    }
    const double rho0 = GetProperties().GetValue(DENSITY);

    // Set the auxiliary references matching the automatic differentiation symbols
    const auto& N = kinematic_variables.N;
    const auto& DN = kinematic_variables.DN_DX;
    const auto& u = kinematic_variables.Displacements;
    const auto& th = kinematic_variables.JacobianDeterminant;
    const auto& S = constitutive_variables.StressVector;
    const auto& C = constitutive_variables.ConstitutiveMatrix;

    // Aux RHS and LHS
    //TODO: To be removed
    BoundedMatrix<double,LocalSize,LocalSize> lhs;

    const SizeType n_gauss = r_geometry.IntegrationPointsNumber(GetIntegrationMethod());
    const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());
    for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
        // Calculate kinematics
        CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());
        const double w_gauss = kinematic_variables.detJ0 * r_integration_points[i_gauss].Weight();

        // Calculate the constitutive response
        CalculateConstitutiveVariables(kinematic_variables, constitutive_variables, cons_law_values, i_gauss, r_geometry.IntegrationPoints(this->GetIntegrationMethod()), ConstitutiveLaw::StressMeasure_Cauchy);

        // Calculate and add the LHS Gauss point contributions
        const double clhs0 = DN(0,1)*u(0,0);
const double clhs1 = DN(1,1)*u(1,0);
const double clhs2 = DN(2,1)*u(2,0);
const double clhs3 = clhs0 + clhs1 + clhs2;
const double clhs4 = DN(0,1)*clhs3;
const double clhs5 = C(0,1)*clhs4;
const double clhs6 = DN(0,0)*u(0,0);
const double clhs7 = DN(1,0)*u(1,0);
const double clhs8 = DN(2,0)*u(2,0);
const double clhs9 = clhs6 + clhs7 + clhs8 + 1;
const double clhs10 = DN(0,0)*clhs9;
const double clhs11 = C(0,0)*clhs10;
const double clhs12 = DN(0,0)*clhs3 + DN(0,1)*clhs9;
const double clhs13 = C(0,2)*clhs12;
const double clhs14 = 1.0*clhs11 + clhs13 + 1.0*clhs5;
const double clhs15 = N[0]*th[0] + N[1]*th[1] + N[2]*th[2];
const double clhs16 = 1.0/clhs15;
const double clhs17 = clhs16*clhs9;
const double clhs18 = C(1,2)*clhs4;
const double clhs19 = C(0,2)*clhs10;
const double clhs20 = C(2,2)*clhs12;
const double clhs21 = 1.0*clhs18 + 1.0*clhs19 + clhs20;
const double clhs22 = clhs16*clhs3;
const double clhs23 = DN(0,0)*u(0,1);
const double clhs24 = DN(1,0)*u(1,1);
const double clhs25 = DN(2,0)*u(2,1);
const double clhs26 = clhs23 + clhs24 + clhs25;
const double clhs27 = pow(clhs26, 2);
const double clhs28 = pow(clhs9, 2);
const double clhs29 = 0.5*clhs16*clhs27 + 0.5*clhs16*clhs28 - 0.5;
const double clhs30 = pow(clhs3, 2);
const double clhs31 = DN(0,1)*u(0,1);
const double clhs32 = DN(1,1)*u(1,1);
const double clhs33 = DN(2,1)*u(2,1);
const double clhs34 = clhs31 + clhs32 + clhs33;
const double clhs35 = clhs34 + 1;
const double clhs36 = pow(clhs35, 2);
const double clhs37 = 0.5*clhs16*clhs30 + 0.5*clhs16*clhs36 - 0.5;
const double clhs38 = clhs26*clhs35 + clhs3*clhs9;
const double clhs39 = clhs16*clhs38;
const double clhs40 = C(0,0)*clhs29 + C(0,1)*clhs37 + C(0,2)*clhs39;
const double clhs41 = C(0,2)*clhs29 + C(1,2)*clhs37 + C(2,2)*clhs39;
const double clhs42 = DN(0,0)*clhs40 + DN(0,1)*clhs41;
const double clhs43 = clhs14*clhs17 + clhs21*clhs22 + clhs42;
const double clhs44 = C(1,1)*clhs4;
const double clhs45 = C(0,1)*clhs10;
const double clhs46 = C(1,2)*clhs12;
const double clhs47 = 1.0*clhs44 + 1.0*clhs45 + clhs46;
const double clhs48 = C(0,1)*clhs29 + C(1,1)*clhs37 + C(1,2)*clhs39;
const double clhs49 = DN(0,0)*clhs41 + DN(0,1)*clhs48;
const double clhs50 = clhs17*clhs21 + clhs22*clhs47 + clhs49;
const double clhs51 = DN(0,0)*clhs26;
const double clhs52 = C(0,2)*clhs51;
const double clhs53 = DN(0,1)*clhs35;
const double clhs54 = C(1,2)*clhs53;
const double clhs55 = DN(0,0)*clhs35 + DN(0,1)*clhs26;
const double clhs56 = C(2,2)*clhs55;
const double clhs57 = 1.0*clhs52 + 1.0*clhs54 + clhs56;
const double clhs58 = C(0,0)*clhs51;
const double clhs59 = C(0,1)*clhs53;
const double clhs60 = C(0,2)*clhs55;
const double clhs61 = 1.0*clhs58 + 1.0*clhs59 + clhs60;
const double clhs62 = clhs3*clhs57 + clhs61*clhs9;
const double clhs63 = C(0,1)*clhs51;
const double clhs64 = C(1,1)*clhs53;
const double clhs65 = C(1,2)*clhs55;
const double clhs66 = 1.0*clhs63 + 1.0*clhs64 + clhs65;
const double clhs67 = clhs3*clhs66 + clhs57*clhs9;
const double clhs68 = pow(clhs15, -2.0);
const double clhs69 = clhs27 + clhs28;
const double clhs70 = C(0,2)*clhs69;
const double clhs71 = clhs30 + clhs36;
const double clhs72 = C(1,2)*clhs71;
const double clhs73 = 1.0*clhs38;
const double clhs74 = C(2,2)*clhs73 + 0.5*clhs70 + 0.5*clhs72;
const double clhs75 = C(0,0)*clhs69;
const double clhs76 = C(0,1)*clhs71;
const double clhs77 = C(0,2)*clhs73 + 0.5*clhs75 + 0.5*clhs76;
const double clhs78 = clhs3*clhs74 + clhs77*clhs9;
const double clhs79 = C(0,1)*clhs69;
const double clhs80 = C(1,1)*clhs71;
const double clhs81 = C(1,2)*clhs73 + 0.5*clhs79 + 0.5*clhs80;
const double clhs82 = clhs3*clhs81 + clhs74*clhs9;
const double clhs83 = clhs68*(DN(0,0)*clhs78 + DN(0,1)*clhs82);
const double clhs84 = DN(1,1)*clhs3;
const double clhs85 = C(0,1)*clhs84;
const double clhs86 = DN(1,0)*clhs9;
const double clhs87 = C(0,0)*clhs86;
const double clhs88 = DN(1,0)*clhs3 + DN(1,1)*clhs9;
const double clhs89 = C(0,2)*clhs88;
const double clhs90 = 1.0*clhs85 + 1.0*clhs87 + clhs89;
const double clhs91 = C(1,2)*clhs84;
const double clhs92 = C(0,2)*clhs86;
const double clhs93 = C(2,2)*clhs88;
const double clhs94 = 1.0*clhs91 + 1.0*clhs92 + clhs93;
const double clhs95 = DN(1,0)*clhs40 + DN(1,1)*clhs41;
const double clhs96 = clhs17*clhs90 + clhs22*clhs94 + clhs95;
const double clhs97 = C(1,1)*clhs84;
const double clhs98 = C(0,1)*clhs86;
const double clhs99 = C(1,2)*clhs88;
const double clhs100 = 1.0*clhs97 + 1.0*clhs98 + clhs99;
const double clhs101 = DN(1,0)*clhs41 + DN(1,1)*clhs48;
const double clhs102 = clhs100*clhs22 + clhs101 + clhs17*clhs94;
const double clhs103 = DN(1,0)*clhs26;
const double clhs104 = C(0,2)*clhs103;
const double clhs105 = DN(1,1)*clhs35;
const double clhs106 = C(1,2)*clhs105;
const double clhs107 = DN(1,0)*clhs35 + DN(1,1)*clhs26;
const double clhs108 = C(2,2)*clhs107;
const double clhs109 = 1.0*clhs104 + 1.0*clhs106 + clhs108;
const double clhs110 = C(0,0)*clhs103;
const double clhs111 = C(0,1)*clhs105;
const double clhs112 = C(0,2)*clhs107;
const double clhs113 = 1.0*clhs110 + 1.0*clhs111 + clhs112;
const double clhs114 = clhs109*clhs3 + clhs113*clhs9;
const double clhs115 = C(0,1)*clhs103;
const double clhs116 = C(1,1)*clhs105;
const double clhs117 = C(1,2)*clhs107;
const double clhs118 = 1.0*clhs115 + 1.0*clhs116 + clhs117;
const double clhs119 = clhs109*clhs9 + clhs118*clhs3;
const double clhs120 = DN(2,1)*clhs3;
const double clhs121 = C(0,1)*clhs120;
const double clhs122 = DN(2,0)*clhs9;
const double clhs123 = C(0,0)*clhs122;
const double clhs124 = DN(2,0)*clhs3 + DN(2,1)*clhs9;
const double clhs125 = C(0,2)*clhs124;
const double clhs126 = 1.0*clhs121 + 1.0*clhs123 + clhs125;
const double clhs127 = C(1,2)*clhs120;
const double clhs128 = C(0,2)*clhs122;
const double clhs129 = C(2,2)*clhs124;
const double clhs130 = 1.0*clhs127 + 1.0*clhs128 + clhs129;
const double clhs131 = DN(2,0)*clhs40 + DN(2,1)*clhs41;
const double clhs132 = clhs126*clhs17 + clhs130*clhs22 + clhs131;
const double clhs133 = C(1,1)*clhs120;
const double clhs134 = C(0,1)*clhs122;
const double clhs135 = C(1,2)*clhs124;
const double clhs136 = 1.0*clhs133 + 1.0*clhs134 + clhs135;
const double clhs137 = DN(2,0)*clhs41 + DN(2,1)*clhs48;
const double clhs138 = clhs130*clhs17 + clhs136*clhs22 + clhs137;
const double clhs139 = DN(2,0)*clhs26;
const double clhs140 = C(0,2)*clhs139;
const double clhs141 = DN(2,1)*clhs35;
const double clhs142 = C(1,2)*clhs141;
const double clhs143 = DN(2,0)*clhs35 + DN(2,1)*clhs26;
const double clhs144 = C(2,2)*clhs143;
const double clhs145 = 1.0*clhs140 + 1.0*clhs142 + clhs144;
const double clhs146 = C(0,0)*clhs139;
const double clhs147 = C(0,1)*clhs141;
const double clhs148 = C(0,2)*clhs143;
const double clhs149 = 1.0*clhs146 + 1.0*clhs147 + clhs148;
const double clhs150 = clhs145*clhs3 + clhs149*clhs9;
const double clhs151 = C(0,1)*clhs139;
const double clhs152 = C(1,1)*clhs141;
const double clhs153 = C(1,2)*clhs143;
const double clhs154 = 1.0*clhs151 + 1.0*clhs152 + clhs153;
const double clhs155 = clhs145*clhs9 + clhs154*clhs3;
const double clhs156 = clhs14*clhs26 + clhs21*clhs35;
const double clhs157 = clhs21*clhs26 + clhs35*clhs47;
const double clhs158 = clhs16*clhs26;
const double clhs159 = clhs16*clhs35;
const double clhs160 = clhs158*clhs61 + clhs159*clhs57 + clhs42;
const double clhs161 = clhs158*clhs57 + clhs159*clhs66 + clhs49;
const double clhs162 = clhs26*clhs77 + clhs35*clhs74;
const double clhs163 = clhs26*clhs74 + clhs35*clhs81;
const double clhs164 = clhs68*(DN(0,0)*clhs162 + DN(0,1)*clhs163);
const double clhs165 = clhs26*clhs90 + clhs35*clhs94;
const double clhs166 = clhs100*clhs35 + clhs26*clhs94;
const double clhs167 = clhs109*clhs159 + clhs113*clhs158 + clhs95;
const double clhs168 = clhs101 + clhs109*clhs158 + clhs118*clhs159;
const double clhs169 = clhs126*clhs26 + clhs130*clhs35;
const double clhs170 = clhs130*clhs26 + clhs136*clhs35;
const double clhs171 = clhs131 + clhs145*clhs159 + clhs149*clhs158;
const double clhs172 = clhs137 + clhs145*clhs158 + clhs154*clhs159;
const double clhs173 = DN(0,0)*clhs32 + DN(0,0)*clhs33 + DN(0,0) - DN(0,1)*clhs24 - DN(0,1)*clhs25;
const double clhs174 = clhs11 + clhs13 + clhs5;
const double clhs175 = clhs18 + clhs19 + clhs20;
const double clhs176 = 2*clhs38;
const double clhs177 = C(0,2)*clhs176 + clhs75 + clhs76;
const double clhs178 = DN(0,0)*clhs177;
const double clhs179 = -clhs0*clhs24 - clhs0*clhs25 - clhs1*clhs23 - clhs1*clhs25 - clhs2*clhs23 - clhs2*clhs24 + clhs31*clhs7 + clhs31*clhs8 + clhs32*clhs6 + clhs32*clhs8 + clhs33*clhs6 + clhs33*clhs7 + clhs34 + clhs9;
const double clhs180 = 0.25/clhs179;
const double clhs181 = clhs173*clhs180;
const double clhs182 = C(2,2)*clhs176 + clhs70 + clhs72;
const double clhs183 = DN(0,1)*clhs182;
const double clhs184 = DN(0,0)*th[0] + DN(1,0)*th[1] + DN(2,0)*th[2];
const double clhs185 = sqrt(clhs179/clhs15);
const double clhs186 = clhs185*tau;
const double clhs187 = clhs16*clhs186;
const double clhs188 = clhs184*clhs187;
const double clhs189 = clhs44 + clhs45 + clhs46;
const double clhs190 = DN(0,0)*clhs182;
const double clhs191 = C(1,2)*clhs176 + clhs79 + clhs80;
const double clhs192 = DN(0,1)*clhs191;
const double clhs193 = DN(0,1)*th[0] + DN(1,1)*th[1] + DN(2,1)*th[2];
const double clhs194 = clhs187*clhs193;
const double clhs195 = -DN(0,0)*clhs1 - DN(0,0)*clhs2 + DN(0,1)*clhs7 + DN(0,1)*clhs8 + DN(0,1);
const double clhs196 = clhs58 + clhs59 + clhs60;
const double clhs197 = clhs52 + clhs54 + clhs56;
const double clhs198 = clhs180*clhs195;
const double clhs199 = clhs63 + clhs64 + clhs65;
const double clhs200 = (1.0/2.0)*clhs187;
const double clhs201 = clhs200*(clhs178 + clhs183);
const double clhs202 = clhs200*(clhs190 + clhs192);
const double clhs203 = (1.0/2.0)*clhs186*clhs68;
const double clhs204 = N[0]*clhs203;
const double clhs205 = 2.0*clhs38;
const double clhs206 = C(0,2)*clhs205 + 1.0*clhs75 + 1.0*clhs76;
const double clhs207 = C(2,2)*clhs205 + 1.0*clhs70 + 1.0*clhs72;
const double clhs208 = clhs184*(DN(0,0)*clhs206 + DN(0,1)*clhs207 + 0.5*clhs178 + 0.5*clhs183);
const double clhs209 = C(1,2)*clhs205 + 1.0*clhs79 + 1.0*clhs80;
const double clhs210 = clhs193*(DN(0,0)*clhs207 + DN(0,1)*clhs209 + 0.5*clhs190 + 0.5*clhs192);
const double clhs211 = N[0]*b(0,1) + N[1]*b(1,1) + N[2]*b(2,1);
const double clhs212 = rho0*tau;
const double clhs213 = clhs212*(DN(0,0)*DN(1,1) - DN(0,1)*DN(1,0));
const double clhs214 = clhs211*clhs213;
const double clhs215 = DN(1,0)*clhs31 + DN(1,0)*clhs33 + DN(1,0) - DN(1,1)*clhs23 - DN(1,1)*clhs25;
const double clhs216 = clhs85 + clhs87 + clhs89;
const double clhs217 = clhs91 + clhs92 + clhs93;
const double clhs218 = clhs180*clhs215;
const double clhs219 = clhs97 + clhs98 + clhs99;
const double clhs220 = N[0]*b(0,0) + N[1]*b(1,0) + N[2]*b(2,0);
const double clhs221 = clhs213*clhs220;
const double clhs222 = -DN(1,0)*clhs0 - DN(1,0)*clhs2 + DN(1,1)*clhs6 + DN(1,1)*clhs8 + DN(1,1);
const double clhs223 = clhs110 + clhs111 + clhs112;
const double clhs224 = clhs104 + clhs106 + clhs108;
const double clhs225 = clhs180*clhs222;
const double clhs226 = clhs115 + clhs116 + clhs117;
const double clhs227 = N[0]*N[1];
const double clhs228 = N[1]*clhs203;
const double clhs229 = clhs212*(DN(0,0)*DN(2,1) - DN(0,1)*DN(2,0));
const double clhs230 = clhs211*clhs229;
const double clhs231 = DN(2,0)*clhs31 + DN(2,0)*clhs32 + DN(2,0) - DN(2,1)*clhs23 - DN(2,1)*clhs24;
const double clhs232 = clhs121 + clhs123 + clhs125;
const double clhs233 = clhs127 + clhs128 + clhs129;
const double clhs234 = clhs180*clhs231;
const double clhs235 = clhs133 + clhs134 + clhs135;
const double clhs236 = clhs220*clhs229;
const double clhs237 = -DN(2,0)*clhs0 - DN(2,0)*clhs1 + DN(2,1)*clhs6 + DN(2,1)*clhs7 + DN(2,1);
const double clhs238 = clhs146 + clhs147 + clhs148;
const double clhs239 = clhs140 + clhs142 + clhs144;
const double clhs240 = clhs180*clhs237;
const double clhs241 = clhs151 + clhs152 + clhs153;
const double clhs242 = N[0]*N[2];
const double clhs243 = N[2]*clhs203;
const double clhs244 = clhs68*(DN(1,0)*clhs78 + DN(1,1)*clhs82);
const double clhs245 = clhs68*(DN(1,0)*clhs162 + DN(1,1)*clhs163);
const double clhs246 = DN(1,0)*clhs177;
const double clhs247 = DN(1,1)*clhs182;
const double clhs248 = DN(1,0)*clhs182;
const double clhs249 = DN(1,1)*clhs191;
const double clhs250 = clhs200*(clhs246 + clhs247);
const double clhs251 = clhs200*(clhs248 + clhs249);
const double clhs252 = clhs184*(DN(1,0)*clhs206 + DN(1,1)*clhs207 + 0.5*clhs246 + 0.5*clhs247);
const double clhs253 = clhs193*(DN(1,0)*clhs207 + DN(1,1)*clhs209 + 0.5*clhs248 + 0.5*clhs249);
const double clhs254 = clhs212*(DN(1,0)*DN(2,1) - DN(1,1)*DN(2,0));
const double clhs255 = clhs211*clhs254;
const double clhs256 = clhs220*clhs254;
const double clhs257 = N[1]*N[2];
const double clhs258 = clhs68*(DN(2,0)*clhs78 + DN(2,1)*clhs82);
const double clhs259 = clhs68*(DN(2,0)*clhs162 + DN(2,1)*clhs163);
const double clhs260 = DN(2,0)*clhs177;
const double clhs261 = DN(2,1)*clhs182;
const double clhs262 = DN(2,0)*clhs182;
const double clhs263 = DN(2,1)*clhs191;
const double clhs264 = clhs200*(clhs260 + clhs261);
const double clhs265 = clhs200*(clhs262 + clhs263);
const double clhs266 = clhs184*(DN(2,0)*clhs206 + DN(2,1)*clhs207 + 0.5*clhs260 + 0.5*clhs261);
const double clhs267 = clhs193*(DN(2,0)*clhs207 + DN(2,1)*clhs209 + 0.5*clhs262 + 0.5*clhs263);
lhs(0,0)=-DN(0,0)*clhs43 - DN(0,1)*clhs50;
lhs(0,1)=-clhs16*(DN(0,0)*clhs62 + DN(0,1)*clhs67);
lhs(0,2)=N[0]*clhs83;
lhs(0,3)=-DN(0,0)*clhs96 - DN(0,1)*clhs102;
lhs(0,4)=-clhs16*(DN(0,0)*clhs114 + DN(0,1)*clhs119);
lhs(0,5)=N[1]*clhs83;
lhs(0,6)=-DN(0,0)*clhs132 - DN(0,1)*clhs138;
lhs(0,7)=-clhs16*(DN(0,0)*clhs150 + DN(0,1)*clhs155);
lhs(0,8)=N[2]*clhs83;
lhs(1,0)=-clhs16*(DN(0,0)*clhs156 + DN(0,1)*clhs157);
lhs(1,1)=-DN(0,0)*clhs160 - DN(0,1)*clhs161;
lhs(1,2)=N[0]*clhs164;
lhs(1,3)=-clhs16*(DN(0,0)*clhs165 + DN(0,1)*clhs166);
lhs(1,4)=-DN(0,0)*clhs167 - DN(0,1)*clhs168;
lhs(1,5)=N[1]*clhs164;
lhs(1,6)=-clhs16*(DN(0,0)*clhs169 + DN(0,1)*clhs170);
lhs(1,7)=-DN(0,0)*clhs171 - DN(0,1)*clhs172;
lhs(1,8)=N[2]*clhs164;
lhs(2,0)=-N[0]*clhs173 + clhs188*(DN(0,0)*clhs174 + DN(0,1)*clhs175 + clhs178*clhs181 + clhs181*clhs183) + clhs194*(DN(0,0)*clhs175 + DN(0,1)*clhs189 + clhs181*clhs190 + clhs181*clhs192);
lhs(2,1)=-N[0]*clhs195 + clhs188*(DN(0,0)*clhs196 + DN(0,1)*clhs197 + clhs178*clhs198 + clhs183*clhs198) + clhs194*(DN(0,0)*clhs197 + DN(0,1)*clhs199 + clhs190*clhs198 + clhs192*clhs198);
lhs(2,2)=DN(0,0)*clhs201 + DN(0,1)*clhs202 + pow(N[0], 2) - clhs204*clhs208 - clhs204*clhs210;
lhs(2,3)=-N[0]*clhs215 + clhs16*clhs184*clhs185*tau*(DN(0,0)*clhs216 + DN(0,1)*clhs217 + clhs178*clhs218 + clhs183*clhs218) + clhs16*clhs185*clhs193*tau*(DN(0,0)*clhs217 + DN(0,1)*clhs219 + clhs190*clhs218 + clhs192*clhs218) - clhs214;
lhs(2,4)=-N[0]*clhs222 + clhs188*(DN(0,0)*clhs223 + DN(0,1)*clhs224 + clhs178*clhs225 + clhs183*clhs225) + clhs194*(DN(0,0)*clhs224 + DN(0,1)*clhs226 + clhs190*clhs225 + clhs192*clhs225) + clhs221;
lhs(2,5)=DN(1,0)*clhs201 + DN(1,1)*clhs202 - clhs208*clhs228 - clhs210*clhs228 + clhs227;
lhs(2,6)=-N[0]*clhs231 + clhs16*clhs184*clhs185*tau*(DN(0,0)*clhs232 + DN(0,1)*clhs233 + clhs178*clhs234 + clhs183*clhs234) + clhs16*clhs185*clhs193*tau*(DN(0,0)*clhs233 + DN(0,1)*clhs235 + clhs190*clhs234 + clhs192*clhs234) - clhs230;
lhs(2,7)=-N[0]*clhs237 + clhs188*(DN(0,0)*clhs238 + DN(0,1)*clhs239 + clhs178*clhs240 + clhs183*clhs240) + clhs194*(DN(0,0)*clhs239 + DN(0,1)*clhs241 + clhs190*clhs240 + clhs192*clhs240) + clhs236;
lhs(2,8)=DN(2,0)*clhs201 + DN(2,1)*clhs202 - clhs208*clhs243 - clhs210*clhs243 + clhs242;
lhs(3,0)=-DN(1,0)*clhs43 - DN(1,1)*clhs50;
lhs(3,1)=-clhs16*(DN(1,0)*clhs62 + DN(1,1)*clhs67);
lhs(3,2)=N[0]*clhs244;
lhs(3,3)=-DN(1,0)*clhs96 - DN(1,1)*clhs102;
lhs(3,4)=-clhs16*(DN(1,0)*clhs114 + DN(1,1)*clhs119);
lhs(3,5)=N[1]*clhs244;
lhs(3,6)=-DN(1,0)*clhs132 - DN(1,1)*clhs138;
lhs(3,7)=-clhs16*(DN(1,0)*clhs150 + DN(1,1)*clhs155);
lhs(3,8)=N[2]*clhs244;
lhs(4,0)=-clhs16*(DN(1,0)*clhs156 + DN(1,1)*clhs157);
lhs(4,1)=-DN(1,0)*clhs160 - DN(1,1)*clhs161;
lhs(4,2)=N[0]*clhs245;
lhs(4,3)=-clhs16*(DN(1,0)*clhs165 + DN(1,1)*clhs166);
lhs(4,4)=-DN(1,0)*clhs167 - DN(1,1)*clhs168;
lhs(4,5)=N[1]*clhs245;
lhs(4,6)=-clhs16*(DN(1,0)*clhs169 + DN(1,1)*clhs170);
lhs(4,7)=-DN(1,0)*clhs171 - DN(1,1)*clhs172;
lhs(4,8)=N[2]*clhs245;
lhs(5,0)=-N[1]*clhs173 + clhs188*(DN(1,0)*clhs174 + DN(1,1)*clhs175 + clhs181*clhs246 + clhs181*clhs247) + clhs194*(DN(1,0)*clhs175 + DN(1,1)*clhs189 + clhs181*clhs248 + clhs181*clhs249) + clhs214;
lhs(5,1)=-N[1]*clhs195 + clhs16*clhs184*clhs185*tau*(DN(1,0)*clhs196 + DN(1,1)*clhs197 + clhs198*clhs246 + clhs198*clhs247) + clhs16*clhs185*clhs193*tau*(DN(1,0)*clhs197 + DN(1,1)*clhs199 + clhs198*clhs248 + clhs198*clhs249) - clhs221;
lhs(5,2)=DN(0,0)*clhs250 + DN(0,1)*clhs251 - clhs204*clhs252 - clhs204*clhs253 + clhs227;
lhs(5,3)=-N[1]*clhs215 + clhs188*(DN(1,0)*clhs216 + DN(1,1)*clhs217 + clhs218*clhs246 + clhs218*clhs247) + clhs194*(DN(1,0)*clhs217 + DN(1,1)*clhs219 + clhs218*clhs248 + clhs218*clhs249);
lhs(5,4)=-N[1]*clhs222 + clhs188*(DN(1,0)*clhs223 + DN(1,1)*clhs224 + clhs225*clhs246 + clhs225*clhs247) + clhs194*(DN(1,0)*clhs224 + DN(1,1)*clhs226 + clhs225*clhs248 + clhs225*clhs249);
lhs(5,5)=DN(1,0)*clhs250 + DN(1,1)*clhs251 + pow(N[1], 2) - clhs228*clhs252 - clhs228*clhs253;
lhs(5,6)=-N[1]*clhs231 + clhs16*clhs184*clhs185*tau*(DN(1,0)*clhs232 + DN(1,1)*clhs233 + clhs234*clhs246 + clhs234*clhs247) + clhs16*clhs185*clhs193*tau*(DN(1,0)*clhs233 + DN(1,1)*clhs235 + clhs234*clhs248 + clhs234*clhs249) - clhs255;
lhs(5,7)=-N[1]*clhs237 + clhs188*(DN(1,0)*clhs238 + DN(1,1)*clhs239 + clhs240*clhs246 + clhs240*clhs247) + clhs194*(DN(1,0)*clhs239 + DN(1,1)*clhs241 + clhs240*clhs248 + clhs240*clhs249) + clhs256;
lhs(5,8)=DN(2,0)*clhs250 + DN(2,1)*clhs251 - clhs243*clhs252 - clhs243*clhs253 + clhs257;
lhs(6,0)=-DN(2,0)*clhs43 - DN(2,1)*clhs50;
lhs(6,1)=-clhs16*(DN(2,0)*clhs62 + DN(2,1)*clhs67);
lhs(6,2)=N[0]*clhs258;
lhs(6,3)=-DN(2,0)*clhs96 - DN(2,1)*clhs102;
lhs(6,4)=-clhs16*(DN(2,0)*clhs114 + DN(2,1)*clhs119);
lhs(6,5)=N[1]*clhs258;
lhs(6,6)=-DN(2,0)*clhs132 - DN(2,1)*clhs138;
lhs(6,7)=-clhs16*(DN(2,0)*clhs150 + DN(2,1)*clhs155);
lhs(6,8)=N[2]*clhs258;
lhs(7,0)=-clhs16*(DN(2,0)*clhs156 + DN(2,1)*clhs157);
lhs(7,1)=-DN(2,0)*clhs160 - DN(2,1)*clhs161;
lhs(7,2)=N[0]*clhs259;
lhs(7,3)=-clhs16*(DN(2,0)*clhs165 + DN(2,1)*clhs166);
lhs(7,4)=-DN(2,0)*clhs167 - DN(2,1)*clhs168;
lhs(7,5)=N[1]*clhs259;
lhs(7,6)=-clhs16*(DN(2,0)*clhs169 + DN(2,1)*clhs170);
lhs(7,7)=-DN(2,0)*clhs171 - DN(2,1)*clhs172;
lhs(7,8)=N[2]*clhs259;
lhs(8,0)=-N[2]*clhs173 + clhs188*(DN(2,0)*clhs174 + DN(2,1)*clhs175 + clhs181*clhs260 + clhs181*clhs261) + clhs194*(DN(2,0)*clhs175 + DN(2,1)*clhs189 + clhs181*clhs262 + clhs181*clhs263) + clhs230;
lhs(8,1)=-N[2]*clhs195 + clhs16*clhs184*clhs185*tau*(DN(2,0)*clhs196 + DN(2,1)*clhs197 + clhs198*clhs260 + clhs198*clhs261) + clhs16*clhs185*clhs193*tau*(DN(2,0)*clhs197 + DN(2,1)*clhs199 + clhs198*clhs262 + clhs198*clhs263) - clhs236;
lhs(8,2)=DN(0,0)*clhs264 + DN(0,1)*clhs265 - clhs204*clhs266 - clhs204*clhs267 + clhs242;
lhs(8,3)=-N[2]*clhs215 + clhs188*(DN(2,0)*clhs216 + DN(2,1)*clhs217 + clhs218*clhs260 + clhs218*clhs261) + clhs194*(DN(2,0)*clhs217 + DN(2,1)*clhs219 + clhs218*clhs262 + clhs218*clhs263) + clhs255;
lhs(8,4)=-N[2]*clhs222 + clhs16*clhs184*clhs185*tau*(DN(2,0)*clhs223 + DN(2,1)*clhs224 + clhs225*clhs260 + clhs225*clhs261) + clhs16*clhs185*clhs193*tau*(DN(2,0)*clhs224 + DN(2,1)*clhs226 + clhs225*clhs262 + clhs225*clhs263) - clhs256;
lhs(8,5)=DN(1,0)*clhs264 + DN(1,1)*clhs265 - clhs228*clhs266 - clhs228*clhs267 + clhs257;
lhs(8,6)=-N[2]*clhs231 + clhs188*(DN(2,0)*clhs232 + DN(2,1)*clhs233 + clhs234*clhs260 + clhs234*clhs261) + clhs194*(DN(2,0)*clhs233 + DN(2,1)*clhs235 + clhs234*clhs262 + clhs234*clhs263);
lhs(8,7)=-N[2]*clhs237 + clhs188*(DN(2,0)*clhs238 + DN(2,1)*clhs239 + clhs240*clhs260 + clhs240*clhs261) + clhs194*(DN(2,0)*clhs239 + DN(2,1)*clhs241 + clhs240*clhs262 + clhs240*clhs263);
lhs(8,8)=DN(2,0)*clhs264 + DN(2,1)*clhs265 + pow(N[2], 2) - clhs243*clhs266 - clhs243*clhs267;

        //TODO: Amend this once the assembly is done in the input arrays
        rLeftHandSideMatrix += w_gauss * lhs;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void TotalLagrangianMixedDetJElement<3>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geometry = GetGeometry();
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType block_size = dim + 1;
    const SizeType matrix_size = block_size * n_nodes;
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();

    // Check LHS size
    if (rLeftHandSideMatrix.size1() != matrix_size || rLeftHandSideMatrix.size2() != matrix_size) {
        rLeftHandSideMatrix.resize(matrix_size, matrix_size, false);
    }

    // Create the kinematics container and fill the nodal data
    KinematicVariables kinematic_variables;
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
        const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        for (IndexType d = 0; d < dim; ++d) {
            kinematic_variables.Displacements(i_node * dim + d) = r_disp[d];
        }
        kinematic_variables.JacobianDeterminant[i_node] = r_geometry[i_node].FastGetSolutionStepValue(DETERMINANT_F);
    }

    // Create the constitutive variables and values containers
    ConstitutiveVariables constitutive_variables;
    ConstitutiveLaw::Parameters cons_law_values(r_geometry, GetProperties(), rCurrentProcessInfo);
    auto& r_cons_law_options = cons_law_values.GetOptions();
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // Calculate the RHS and LHS contributions
    rLeftHandSideMatrix.clear();

    // Calculate stabilization constant
    const double c_tau = 2.0;
    const double h = ElementSizeCalculator<3,NumNodes>::AverageElementSize(r_geometry);
    const double mu = 1.0; //FIXME: This is the Lame constant. Compute it.
    const double tau = c_tau * std::pow(h,2) / (2.0 * mu);

    // Set data for the body force calculation
    BoundedMatrix<double, NumNodes, 3> b;
    for (IndexType i_node = 0; i_node < NumNodes; ++i_node) {
        const array_1d<double,3>& r_b_i = r_geometry[i_node].FastGetSolutionStepValue(BODY_FORCE);
        for (IndexType d = 0; d < 3; ++d) {
            b(i_node, d) = r_b_i[d];
        }
    }
    const double rho0 = GetProperties().GetValue(DENSITY);

    // Set the auxiliary references matching the automatic differentiation symbols
    const auto& N = kinematic_variables.N;
    const auto& DN = kinematic_variables.DN_DX;
    const auto& u = kinematic_variables.Displacements;
    const auto& th = kinematic_variables.JacobianDeterminant;
    const auto& S = constitutive_variables.StressVector;
    const auto& C = constitutive_variables.ConstitutiveMatrix;

    // Aux RHS and LHS
    //TODO: To be removed
    BoundedMatrix<double,LocalSize,LocalSize> lhs;

    const SizeType n_gauss = r_geometry.IntegrationPointsNumber(GetIntegrationMethod());
    const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());
    for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
        // Calculate kinematics
        CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());
        const double w_gauss = kinematic_variables.detJ0 * r_integration_points[i_gauss].Weight();

        // Calculate the constitutive response
        CalculateConstitutiveVariables(kinematic_variables, constitutive_variables, cons_law_values, i_gauss, r_geometry.IntegrationPoints(this->GetIntegrationMethod()), ConstitutiveLaw::StressMeasure_Cauchy);

        // Calculate and add the LHS Gauss point contributions
        //substitute_lhs_3D_4N
        //TODO: Amend this once the assembly is done in the input arrays
        rLeftHandSideMatrix += w_gauss * lhs;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void TotalLagrangianMixedDetJElement<2>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geometry = GetGeometry();
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType block_size = dim + 1;
    const SizeType matrix_size = block_size * n_nodes;
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();

    // Check RHS size
    if (rRightHandSideVector.size() != matrix_size) {
        rRightHandSideVector.resize(matrix_size, false);
    }

    // Create the kinematics container and fill the nodal data
    KinematicVariables kinematic_variables;
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
        const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        for (IndexType d = 0; d < dim; ++d) {
            kinematic_variables.Displacements(i_node * dim + d) = r_disp[d];
        }
        kinematic_variables.JacobianDeterminant[i_node] = r_geometry[i_node].FastGetSolutionStepValue(DETERMINANT_F);
    }

    // Create the constitutive variables and values containers
    ConstitutiveVariables constitutive_variables;
    ConstitutiveLaw::Parameters cons_law_values(r_geometry, GetProperties(), rCurrentProcessInfo);
    auto& r_cons_law_options = cons_law_values.GetOptions();
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // Calculate the RHS and LHS contributions
    rRightHandSideVector.clear();

    // Calculate stabilization constant
    const double c_tau = 2.0;
    const double h = ElementSizeCalculator<2,NumNodes>::AverageElementSize(r_geometry);
    const double mu = 1.0; //FIXME: This is the Lame constant. Compute it.
    const double tau = c_tau * std::pow(h,2) / (2.0 * mu);

    // Set data for the body force calculation
    BoundedMatrix<double, NumNodes, 2> b;
    for (IndexType i_node = 0; i_node < NumNodes; ++i_node) {
        const array_1d<double,3>& r_b_i = r_geometry[i_node].FastGetSolutionStepValue(BODY_FORCE);
        for (IndexType d = 0; d < 2; ++d) {
            b(i_node, d) = r_b_i[d];
        }
    }
    const double rho0 = GetProperties().GetValue(DENSITY);

    // Set the auxiliary references matching the automatic differentiation symbols
    const auto& N = kinematic_variables.N;
    const auto& DN = kinematic_variables.DN_DX;
    const auto& u = kinematic_variables.Displacements;
    const auto& th = kinematic_variables.JacobianDeterminant;
    const auto& S = constitutive_variables.StressVector;
    const auto& C = constitutive_variables.ConstitutiveMatrix;

    // Aux RHS and LHS
    //TODO: To be removed
    BoundedVector<double,LocalSize> rhs;

    const SizeType n_gauss = r_geometry.IntegrationPointsNumber(GetIntegrationMethod());
    const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());
    for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
        // Calculate kinematics
        CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());
        const double w_gauss = kinematic_variables.detJ0 * r_integration_points[i_gauss].Weight();

        // Calculate the constitutive response
        CalculateConstitutiveVariables(kinematic_variables, constitutive_variables, cons_law_values, i_gauss, r_geometry.IntegrationPoints(this->GetIntegrationMethod()), ConstitutiveLaw::StressMeasure_Cauchy);

        // Calculate and add the RHS Gauss point contribution
        const double crhs0 = N[0]*b(0,0) + N[1]*b(1,0) + N[2]*b(2,0);
const double crhs1 = N[0]*rho0;
const double crhs2 = DN(0,1)*u(0,0);
const double crhs3 = DN(1,1)*u(1,0);
const double crhs4 = DN(2,1)*u(2,0);
const double crhs5 = crhs2 + crhs3 + crhs4;
const double crhs6 = DN(0,0)*u(0,0);
const double crhs7 = DN(1,0)*u(1,0);
const double crhs8 = DN(2,0)*u(2,0);
const double crhs9 = crhs6 + crhs7 + crhs8 + 1;
const double crhs10 = S[0]*crhs9 + S[2]*crhs5;
const double crhs11 = S[1]*crhs5 + S[2]*crhs9;
const double crhs12 = N[0]*b(0,1) + N[1]*b(1,1) + N[2]*b(2,1);
const double crhs13 = DN(0,0)*u(0,1);
const double crhs14 = DN(1,0)*u(1,1);
const double crhs15 = DN(2,0)*u(2,1);
const double crhs16 = crhs13 + crhs14 + crhs15;
const double crhs17 = DN(0,1)*u(0,1);
const double crhs18 = DN(1,1)*u(1,1);
const double crhs19 = DN(2,1)*u(2,1);
const double crhs20 = crhs17 + crhs18 + crhs19;
const double crhs21 = crhs20 + 1;
const double crhs22 = S[0]*crhs16 + S[2]*crhs21;
const double crhs23 = S[1]*crhs21 + S[2]*crhs16;
const double crhs24 = rho0*tau;
const double crhs25 = crhs12*crhs24;
const double crhs26 = crhs0*crhs24;
const double crhs27 = N[0]*th[0];
const double crhs28 = N[1]*th[1];
const double crhs29 = N[2]*th[2];
const double crhs30 = -crhs13*crhs3 - crhs13*crhs4 - crhs14*crhs2 - crhs14*crhs4 - crhs15*crhs2 - crhs15*crhs3 + crhs17*crhs7 + crhs17*crhs8 + crhs18*crhs6 + crhs18*crhs8 + crhs19*crhs6 + crhs19*crhs7 + crhs20 + crhs9;
const double crhs31 = -crhs27 - crhs28 - crhs29 + crhs30;
const double crhs32 = pow(crhs16, 2) + pow(crhs9, 2);
const double crhs33 = pow(crhs21, 2) + pow(crhs5, 2);
const double crhs34 = 2*crhs16*crhs21 + 2*crhs5*crhs9;
const double crhs35 = C(0,0)*crhs32 + C(0,1)*crhs33 + C(0,2)*crhs34;
const double crhs36 = C(0,2)*crhs32 + C(1,2)*crhs33 + C(2,2)*crhs34;
const double crhs37 = crhs27 + crhs28 + crhs29;
const double crhs38 = (1.0/2.0)*1.0/crhs37*tau*sqrt(crhs30/crhs37);
const double crhs39 = crhs38*(DN(0,0)*th[0] + DN(1,0)*th[1] + DN(2,0)*th[2]);
const double crhs40 = C(0,1)*crhs32 + C(1,1)*crhs33 + C(1,2)*crhs34;
const double crhs41 = crhs38*(DN(0,1)*th[0] + DN(1,1)*th[1] + DN(2,1)*th[2]);
const double crhs42 = N[1]*rho0;
const double crhs43 = N[2]*rho0;
rhs[0]=DN(0,0)*crhs10 + DN(0,1)*crhs11 - crhs0*crhs1;
rhs[1]=DN(0,0)*crhs22 + DN(0,1)*crhs23 - crhs1*crhs12;
rhs[2]=N[0]*crhs31 - crhs25*(-DN(0,0)*crhs5 + DN(0,1)*crhs9) - crhs26*(DN(0,0)*crhs21 - DN(0,1)*crhs16) - crhs39*(DN(0,0)*crhs35 + DN(0,1)*crhs36) - crhs41*(DN(0,0)*crhs36 + DN(0,1)*crhs40);
rhs[3]=DN(1,0)*crhs10 + DN(1,1)*crhs11 - crhs0*crhs42;
rhs[4]=DN(1,0)*crhs22 + DN(1,1)*crhs23 - crhs12*crhs42;
rhs[5]=N[1]*crhs31 - crhs25*(-DN(1,0)*crhs5 + DN(1,1)*crhs9) - crhs26*(DN(1,0)*crhs21 - DN(1,1)*crhs16) - crhs39*(DN(1,0)*crhs35 + DN(1,1)*crhs36) - crhs41*(DN(1,0)*crhs36 + DN(1,1)*crhs40);
rhs[6]=DN(2,0)*crhs10 + DN(2,1)*crhs11 - crhs0*crhs43;
rhs[7]=DN(2,0)*crhs22 + DN(2,1)*crhs23 - crhs12*crhs43;
rhs[8]=N[2]*crhs31 - crhs25*(-DN(2,0)*crhs5 + DN(2,1)*crhs9) - crhs26*(DN(2,0)*crhs21 - DN(2,1)*crhs16) - crhs39*(DN(2,0)*crhs35 + DN(2,1)*crhs36) - crhs41*(DN(2,0)*crhs36 + DN(2,1)*crhs40);

        //TODO: Amend this once the assembly is done in the input arrays
        rRightHandSideVector += w_gauss * rhs;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void TotalLagrangianMixedDetJElement<3>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geometry = GetGeometry();
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType block_size = dim + 1;
    const SizeType matrix_size = block_size * n_nodes;
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();

    // Check RHS size
    if (rRightHandSideVector.size() != matrix_size) {
        rRightHandSideVector.resize(matrix_size, false);
    }

    // Create the kinematics container and fill the nodal data
    KinematicVariables kinematic_variables;
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
        const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        for (IndexType d = 0; d < dim; ++d) {
            kinematic_variables.Displacements(i_node * dim + d) = r_disp[d];
        }
        kinematic_variables.JacobianDeterminant[i_node] = r_geometry[i_node].FastGetSolutionStepValue(DETERMINANT_F);
    }

    // Create the constitutive variables and values containers
    ConstitutiveVariables constitutive_variables;
    ConstitutiveLaw::Parameters cons_law_values(r_geometry, GetProperties(), rCurrentProcessInfo);
    auto& r_cons_law_options = cons_law_values.GetOptions();
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // Calculate the RHS and LHS contributions
    rRightHandSideVector.clear();

    // Calculate stabilization constant
    const double c_tau = 2.0;
    const double h = ElementSizeCalculator<3,NumNodes>::AverageElementSize(r_geometry);
    const double mu = 1.0; //FIXME: This is the Lame constant. Compute it.
    const double tau = c_tau * std::pow(h,2) / (2.0 * mu);

    // Set data for the body force calculation
    BoundedMatrix<double, NumNodes, 3> b;
    for (IndexType i_node = 0; i_node < NumNodes; ++i_node) {
        const array_1d<double,3>& r_b_i = r_geometry[i_node].FastGetSolutionStepValue(BODY_FORCE);
        for (IndexType d = 0; d < 3; ++d) {
            b(i_node, d) = r_b_i[d];
        }
    }
    const double rho0 = GetProperties().GetValue(DENSITY);

    // Set the auxiliary references matching the automatic differentiation symbols
    const auto& N = kinematic_variables.N;
    const auto& DN = kinematic_variables.DN_DX;
    const auto& u = kinematic_variables.Displacements;
    const auto& th = kinematic_variables.JacobianDeterminant;
    const auto& S = constitutive_variables.StressVector;
    const auto& C = constitutive_variables.ConstitutiveMatrix;

    // Aux RHS and LHS
    //TODO: To be removed
    BoundedVector<double,LocalSize> rhs;
    BoundedMatrix<double,LocalSize,LocalSize> lhs;

    const SizeType n_gauss = r_geometry.IntegrationPointsNumber(GetIntegrationMethod());
    const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());
    for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
        // Calculate kinematics
        CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());
        const double w_gauss = kinematic_variables.detJ0 * r_integration_points[i_gauss].Weight();

        // Calculate the constitutive response
        CalculateConstitutiveVariables(kinematic_variables, constitutive_variables, cons_law_values, i_gauss, r_geometry.IntegrationPoints(this->GetIntegrationMethod()), ConstitutiveLaw::StressMeasure_Cauchy);

        // Calculate and add the RHS Gauss point contribution
        //substitute_rhs_3D_4N
        //TODO: Amend this once the assembly is done in the input arrays
        rRightHandSideVector += w_gauss * rhs;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void TotalLagrangianMixedDetJElement<TDim>::InitializeMaterial()
{
    KRATOS_TRY

    const auto& r_properties = GetProperties();
    if (r_properties[CONSTITUTIVE_LAW] != nullptr) {
        const auto& r_geometry = GetGeometry();
        const auto& r_N_values = r_geometry.ShapeFunctionsValues(mThisIntegrationMethod);
        IndexType aux = 0;
        for (auto &it_gauss_pt : mConstitutiveLawVector) {
            it_gauss_pt = (r_properties[CONSTITUTIVE_LAW])->Clone();
            (it_gauss_pt)->InitializeMaterial(r_properties, r_geometry, row(r_N_values, aux));
            aux++;
        }
    } else {
        KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
bool TotalLagrangianMixedDetJElement<TDim>::UseElementProvidedStrain() const
{
    return true;
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void TotalLagrangianMixedDetJElement<TDim>::SetConstitutiveVariables(
    KinematicVariables& rThisKinematicVariables,
    ConstitutiveVariables& rThisConstitutiveVariables,
    ConstitutiveLaw::Parameters& rValues,
    const IndexType PointNumber,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints) const
{
    // Here we essentially set the input parameters
    rValues.SetShapeFunctionsValues(rThisKinematicVariables.N); // shape functions
    rValues.SetStrainVector(rThisKinematicVariables.EquivalentStrain); // equivalent total strain
    //TODO: Check if these are really required. I think they shouldn't as we're computing the strain in the element
    // rValues.SetDeterminantF(rThisKinematicVariables.detF); // assuming that det(F) is computed somewhere else
    // rValues.SetDeformationGradientF(rThisKinematicVariables.F); // assuming that F is computed somewhere else

    // Here we set the space on which the results shall be written
    rValues.SetStressVector(rThisConstitutiveVariables.StressVector); //F computed somewhere else
    rValues.SetConstitutiveMatrix(rThisConstitutiveVariables.ConstitutiveMatrix); //assuming the determinant is computed somewhere else
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void TotalLagrangianMixedDetJElement<TDim>::CalculateConstitutiveVariables(
    KinematicVariables& rThisKinematicVariables,
    ConstitutiveVariables& rThisConstitutiveVariables,
    ConstitutiveLaw::Parameters& rValues,
    const IndexType PointNumber,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
    const ConstitutiveLaw::StressMeasure ThisStressMeasure) const
{
    // Set the constitutive variables
    SetConstitutiveVariables(rThisKinematicVariables, rThisConstitutiveVariables, rValues, PointNumber, IntegrationPoints);

    // Actually do the computations in the ConstitutiveLaw
    mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(rValues, ThisStressMeasure); //here the calculations are actually done
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void TotalLagrangianMixedDetJElement<TDim>::CalculateKinematicVariables(
    KinematicVariables& rThisKinematicVariables,
    const IndexType PointNumber,
    const GeometryType::IntegrationMethod& rIntegrationMethod) const
{
    const auto& r_geometry = GetGeometry();
    const auto& r_integration_points = r_geometry.IntegrationPoints(rIntegrationMethod);

    // Shape functions
    rThisKinematicVariables.N = r_geometry.ShapeFunctionsValues(rThisKinematicVariables.N, r_integration_points[PointNumber].Coordinates());

    // Calculate the inverse Jacobian
    GeometryUtils::JacobianOnInitialConfiguration(
        r_geometry,
        r_integration_points[PointNumber],
        rThisKinematicVariables.J0);
    MathUtils<double>::InvertMatrix(
        rThisKinematicVariables.J0,
        rThisKinematicVariables.InvJ0,
        rThisKinematicVariables.detJ0);
    KRATOS_ERROR_IF(rThisKinematicVariables.detJ0 < 0.0)
        << "Element ID: " << this->Id() << " is inverted. det(J0) = " << rThisKinematicVariables.detJ0 << std::endl;

    // Calculate the shape functions gradients
    GeometryUtils::ShapeFunctionsGradients(
        r_geometry.ShapeFunctionsLocalGradients(rIntegrationMethod)[PointNumber],
        rThisKinematicVariables.InvJ0,
        rThisKinematicVariables.DN_DX);

    // Calculate the equivalent total strain
    CalculateEquivalentStrain(rThisKinematicVariables);

    // Compute equivalent F
    //TODO: Check if these are really required. I think they shouldn't as we're computing the strain in the element
    // ComputeEquivalentF(rThisKinematicVariables.F, rThisKinematicVariables.EquivalentStrain);
    // rThisKinematicVariables.detF = MathUtils<double>::Det(rThisKinematicVariables.F);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void TotalLagrangianMixedDetJElement<2>::CalculateEquivalentStrain(KinematicVariables& rThisKinematicVariables) const
{
    // Define references to the auxiliary symbols
    const auto& N = rThisKinematicVariables.N;
    const auto& DN = rThisKinematicVariables.DN_DX;
    const auto& u = rThisKinematicVariables.Displacements;
    const auto& th = rThisKinematicVariables.JacobianDeterminant;
    auto& eq_green_strain = rThisKinematicVariables.EquivalentStrain;

    // Fill the equivalent Green strain values
    const double ceq_green_strain0 = 1.0/(N[0]*th[0] + N[1]*th[1] + N[2]*th[2]);
const double ceq_green_strain1 = DN(0,0)*u(0,1) + DN(1,0)*u(1,1) + DN(2,0)*u(2,1);
const double ceq_green_strain2 = DN(0,0)*u(0,0) + DN(1,0)*u(1,0) + DN(2,0)*u(2,0) + 1;
const double ceq_green_strain3 = DN(0,1)*u(0,0) + DN(1,1)*u(1,0) + DN(2,1)*u(2,0);
const double ceq_green_strain4 = DN(0,1)*u(0,1) + DN(1,1)*u(1,1) + DN(2,1)*u(2,1) + 1;
eq_green_strain[0]=0.5*ceq_green_strain0*pow(ceq_green_strain1, 2) + 0.5*ceq_green_strain0*pow(ceq_green_strain2, 2) - 0.5;
eq_green_strain[1]=0.5*ceq_green_strain0*pow(ceq_green_strain3, 2) + 0.5*ceq_green_strain0*pow(ceq_green_strain4, 2) - 0.5;
eq_green_strain[2]=1.0*ceq_green_strain0*(ceq_green_strain1*ceq_green_strain4 + ceq_green_strain2*ceq_green_strain3);

}

/***********************************************************************************/
/***********************************************************************************/

template<>
void TotalLagrangianMixedDetJElement<3>::CalculateEquivalentStrain(KinematicVariables& rThisKinematicVariables) const
{
    // Define references to the auxiliary symbols
    const auto& N = rThisKinematicVariables.N;
    const auto& DN = rThisKinematicVariables.DN_DX;
    const auto& u = rThisKinematicVariables.Displacements;
    const auto& th = rThisKinematicVariables.JacobianDeterminant;
    auto& eq_green_strain = rThisKinematicVariables.EquivalentStrain;

    // Fill the equivalent Green strain values
    //substitute_green_strain_3D_4N
}

// /***********************************************************************************/
// /***********************************************************************************/

// template<std::size_t TDim>
// void TotalLagrangianMixedDetJElement<TDim>::ComputeEquivalentF(
//     Matrix& rF,
//     const Vector& rStrainTensor) const
// {
//     StructuralMechanicsElementUtilities::ComputeEquivalentF(*this, rStrainTensor, rF);
// }

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
int  TotalLagrangianMixedDetJElement<TDim>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    int check = TotalLagrangianMixedDetJElement::BaseType::Check(rCurrentProcessInfo);

    // Base check
    check = StructuralMechanicsElementUtilities::SolidElementCheck(*this, rCurrentProcessInfo, mConstitutiveLawVector);

    // Checking density
    KRATOS_ERROR_IF_NOT(GetProperties().Has(DENSITY)) << "DENSITY has to be provided for the calculation of body force." << std::endl;

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    const auto& r_geometry = this->GetGeometry();
    for ( IndexType i = 0; i < r_geometry.size(); i++ ) {
        const NodeType& r_node = r_geometry[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DETERMINANT_F,r_node)
        KRATOS_CHECK_DOF_IN_NODE(DETERMINANT_F, r_node)
    }

    return check;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void TotalLagrangianMixedDetJElement<TDim>::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const auto& r_geometry = GetGeometry();
    const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());

    const SizeType n_gauss = r_integration_points.size();
    if (rOutput.size() != n_gauss) {
        rOutput.resize(n_gauss);
    }

    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        GetValueOnConstitutiveLaw(rVariable, rOutput);
    } else {
        CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void TotalLagrangianMixedDetJElement<TDim>::CalculateOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const auto& r_geometry = GetGeometry();
    const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());

    const SizeType n_gauss = r_integration_points.size();
    if (rOutput.size() != n_gauss) {
        rOutput.resize(n_gauss);
    }

    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        GetValueOnConstitutiveLaw(rVariable, rOutput);
    } else if (rVariable == CAUCHY_STRESS_VECTOR || rVariable == PK2_STRESS_VECTOR) {
        // Create and initialize element variables:
        const SizeType n_nodes = r_geometry.PointsNumber();
        const SizeType dim = r_geometry.WorkingSpaceDimension();
        const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

        // Create the kinematics container and fill the nodal data
        KinematicVariables kinematic_variables;
        for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
            const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
            for (IndexType d = 0; d < dim; ++d) {
                kinematic_variables.Displacements(i_node, d) = r_disp[d];
            }
            kinematic_variables.JacobianDeterminant[i_node] = r_geometry[i_node].FastGetSolutionStepValue(DETERMINANT_F);
        }

        // Create the constitutive variables and values containers
        ConstitutiveVariables constitutive_variables;
        ConstitutiveLaw::Parameters cons_law_values(r_geometry, GetProperties(), rCurrentProcessInfo);
        auto& r_cons_law_options = cons_law_values.GetOptions();
        r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);


        for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
            // Calculate kinematics
            CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());

            // Call the constitutive law to update material variables
            if( rVariable == CAUCHY_STRESS_VECTOR) {
                // Compute material reponse
                CalculateConstitutiveVariables(kinematic_variables, constitutive_variables, cons_law_values, i_gauss, r_integration_points, ConstitutiveLaw::StressMeasure_Cauchy);
            } else {
                // Compute material reponse
                CalculateConstitutiveVariables(kinematic_variables, constitutive_variables, cons_law_values, i_gauss, r_integration_points, ConstitutiveLaw::StressMeasure_PK2);
            }

            // Check sizes and save the output stress
            if (rOutput[i_gauss].size() != strain_size) {
                rOutput[i_gauss].resize(strain_size, false);
            }
            rOutput[i_gauss] = constitutive_variables.StressVector;
        }
    } else {
        CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
const Parameters TotalLagrangianMixedDetJElement<TDim>::GetSpecifications() const
{
    const Parameters specifications = Parameters(R"({
        "time_integration"           : ["static"],
        "framework"                  : "lagrangian",
        "symmetric_lhs"              : true,
        "positive_definite_lhs"      : true,
        "output"                     : {
            "gauss_point"            : ["CAUCHY_STRESS_VECTOR"],
            "nodal_historical"       : ["DISPLACEMENT","DETERMINANT_F"],
            "nodal_non_historical"   : [],
            "entity"                 : []
        },
        "required_variables"         : ["DISPLACEMENT","DETERMINANT_F"],
        "required_dofs"              : [],
        "flags_used"                 : [],
        "compatible_geometries"      : ["Triangle2D3", "Tetrahedra3D4"],
        "element_integrates_in_time" : true,
        "compatible_constitutive_laws": {
            "type"        : ["PlaneStrain","PlaneStress","ThreeDimensional"],
            "dimension"   : ["2D","3D"],
            "strain_size" : [3,6]
        },
        "required_polynomial_degree_of_geometry" : 1,
        "documentation"   :
            "This element implements a mixed displacement - volumetric strain formulation with Variational MultiScales (VMS) stabilization. This formulation is capable to deal with materials in the incompressible limit as well as with anisotropy."
    })");

    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    if (TDim == 2) {
        std::vector<std::string> dofs_2d({"DISPLACEMENT_X","DISPLACEMENT_Y","DETERMINANT_F"});
        specifications["required_dofs"].SetStringArray(dofs_2d);
    } else {
        std::vector<std::string> dofs_3d({"DISPLACEMENT_X","DISPLACEMENT_Y","DISPLACEMENT_Z","DETERMINANT_F"});
        specifications["required_dofs"].SetStringArray(dofs_3d);
    }

    return specifications;
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void TotalLagrangianMixedDetJElement<TDim>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, TotalLagrangianMixedDetJElement::BaseType);
    int IntMethod = int(this->GetIntegrationMethod());
    rSerializer.save("IntegrationMethod",IntMethod);
    rSerializer.save("ConstitutiveLawVector", mConstitutiveLawVector);
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void TotalLagrangianMixedDetJElement<TDim>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, TotalLagrangianMixedDetJElement::BaseType);
    int IntMethod;
    rSerializer.load("IntegrationMethod",IntMethod);
    mThisIntegrationMethod = IntegrationMethod(IntMethod);
    rSerializer.load("ConstitutiveLawVector", mConstitutiveLawVector);
}

/***********************************************************************************/
/***********************************************************************************/

// Explicit template instantiations
template class TotalLagrangianMixedDetJElement<2>;
template class TotalLagrangianMixedDetJElement<3>;

} // Namespace Kratos
