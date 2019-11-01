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
//

// System includes

// External includes


// Project includes
#include "includes/define.h"
#include "custom_elements/small_displacement_mixed_volumetric_strain_element.h"
#include "utilities/geometry_utilities.h"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

Element::Pointer SmallDisplacementMixedVolumetricStrainElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<SmallDisplacementMixedVolumetricStrainElement>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer SmallDisplacementMixedVolumetricStrainElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<SmallDisplacementMixedVolumetricStrainElement>( NewId, pGeom, pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer SmallDisplacementMixedVolumetricStrainElement::Clone (
    IndexType NewId,
    NodesArrayType const& rThisNodes) const
{
    KRATOS_TRY

    SmallDisplacementMixedVolumetricStrainElement::Pointer p_new_elem = Kratos::make_intrusive<SmallDisplacementMixedVolumetricStrainElement>(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
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

void SmallDisplacementMixedVolumetricStrainElement::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto &r_geometry = GetGeometry();
    const unsigned int n_nodes = r_geometry.PointsNumber();
    const unsigned int dim = r_geometry.WorkingSpaceDimension();
    const unsigned int dof_size = n_nodes*(dim+1);

    if (rResult.size() != dof_size){
        rResult.resize(dof_size);
    }

    const unsigned int disp_pos = r_geometry[0].GetDofPosition(DISPLACEMENT_X);
    const unsigned int eps_vol_pos = r_geometry[0].GetDofPosition(VOLUMETRIC_STRAIN);

    unsigned int aux_index = 0;
    if (dim == 2) {
        for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
            rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_X, disp_pos).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_Y, disp_pos + 1).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(VOLUMETRIC_STRAIN, eps_vol_pos).EquationId();
        }
    } else {
        for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
            KRATOS_WATCH(aux_index)
            rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_X, disp_pos).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_Y, disp_pos + 1).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_Z, disp_pos + 2).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(VOLUMETRIC_STRAIN, eps_vol_pos).EquationId();
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::GetDofList(
    DofsVectorType& rElementalDofList,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto &r_geometry = GetGeometry();
    const unsigned int n_nodes = r_geometry.PointsNumber();
    const unsigned int dim = r_geometry.WorkingSpaceDimension();
    const unsigned int dof_size  = n_nodes*(dim+1);

    if (rElementalDofList.size() != dof_size){
        rElementalDofList.resize(dof_size);
    }

    if (dim == 2) {
        for(unsigned int i = 0; i < n_nodes; ++i) {
            rElementalDofList[i * (dim + 1)] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_X);
            rElementalDofList[i * (dim + 1) + 1] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_Y);
            rElementalDofList[i * (dim + 1) + 2] = this->GetGeometry()[i].pGetDof(VOLUMETRIC_STRAIN);
        }
    } else if (dim == 3) {
        for(unsigned int i = 0; i < n_nodes; ++i){
            rElementalDofList[i * (dim + 1)] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_X);
            rElementalDofList[i * (dim + 1) + 1] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_Y);
            rElementalDofList[i * (dim + 1) + 2] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_Z);
            rElementalDofList[i * (dim + 1) + 3] = this->GetGeometry()[i].pGetDof(VOLUMETRIC_STRAIN);
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::Initialize()
{
    KRATOS_TRY

    // Integration method initialization
    mThisIntegrationMethod = GeometryData::GI_GAUSS_2;
    const auto &r_integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

    // Constitutive Law Vector initialisation
    if (mConstitutiveLawVector.size() != r_integration_points.size()) {
        mConstitutiveLawVector.resize(r_integration_points.size());
    }

    // Initialize material
    InitializeMaterial();

    KRATOS_CATCH( "" )
}

void SmallDisplacementMixedVolumetricStrainElement::InitializeSolutionStep(ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY

    // Set te constitutive law values
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();
    Vector strain(strain_size);
    Vector stress(strain_size);
    Matrix cons_matrix(strain_size, strain_size);
    ConstitutiveLaw::Parameters cons_law_values(GetGeometry(), GetProperties(), rCurrentProcessInfo);
    auto &r_cons_law_options = cons_law_values.GetOptions();
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    cons_law_values.SetStrainVector(strain);
    cons_law_values.SetStressVector(stress);
    cons_law_values.SetConstitutiveMatrix(cons_matrix);

    // Call the initialize material response
    for (IndexType i_gauss = 0; i_gauss < mConstitutiveLawVector.size(); ++i_gauss) {
        // Call the constitutive law to update material variables
        mConstitutiveLawVector[i_gauss]->InitializeMaterialResponseCauchy(cons_law_values);
    }

    KRATOS_CATCH( "" )
}

void SmallDisplacementMixedVolumetricStrainElement::FinalizeSolutionStep(ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY

    // Set te constitutive law values
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();
    Vector strain(strain_size);
    Vector stress(strain_size);
    Matrix cons_matrix(strain_size, strain_size);
    ConstitutiveLaw::Parameters cons_law_values(GetGeometry(), GetProperties(), rCurrentProcessInfo);
    auto &r_cons_law_options = cons_law_values.GetOptions();
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    cons_law_values.SetStrainVector(strain);
    cons_law_values.SetStressVector(stress);
    cons_law_values.SetConstitutiveMatrix(cons_matrix);

    // Call the initialize material response
    for (IndexType i_gauss = 0; i_gauss < mConstitutiveLawVector.size(); ++i_gauss) {
        // Call the constitutive law to update material variables
        mConstitutiveLawVector[i_gauss]->FinalizeMaterialResponseCauchy(cons_law_values);
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    // TODO: IMPLEMENT IN A MORE EFFICIENT MANNER
    CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
    CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    ProcessInfo& rCurrentProcessInfo)
{
    const auto &r_geometry = GetGeometry();
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType block_size = dim + 1;
    const SizeType matrix_size = block_size * n_nodes;
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();

    // Check RHS size
    if (rLeftHandSideMatrix.size1() != matrix_size || rLeftHandSideMatrix.size2() != matrix_size) {
        rLeftHandSideMatrix.resize(matrix_size, matrix_size, false);
    }

    // Create the kinematics container and fill the nodal data
    KinematicVariables kinematic_variables(strain_size, dim, n_nodes);
    for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
        const auto &r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        for (unsigned int d = 0; d < dim; ++d) {
            kinematic_variables.Displacements(i_node * dim + d) = r_disp[d];
        }
        kinematic_variables.VolumetricNodalStrains[i_node] = r_geometry[i_node].FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
    }

    // Create the constitutive variables and values containers
    ConstitutiveVariables constitutive_variables(strain_size);
    ConstitutiveLaw::Parameters cons_law_values(r_geometry, GetProperties(), rCurrentProcessInfo);
    auto &r_cons_law_options = cons_law_values.GetOptions();
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // Calculate the LHS contributions
    rLeftHandSideMatrix.clear();

    const SizeType n_gauss = r_geometry.IntegrationPointsNumber(GetIntegrationMethod());
    const auto &r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());
    for (unsigned int i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
        // Calculate kinematics
        CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());
        const double w_gauss = kinematic_variables.detJ0 * r_integration_points[i_gauss].Weight();

        // Calculate the constitutive response
        CalculateConstitutiveVariables(kinematic_variables, constitutive_variables, cons_law_values, i_gauss, r_geometry.IntegrationPoints(this->GetIntegrationMethod()), ConstitutiveLaw::StressMeasure_Cauchy);

        // Add momentum volume stress contribution
        double vol_strain = 0.0;
        double vol_stress = 0.0;
        for (unsigned int d = 0; d < dim; ++d) {
            vol_strain += cons_law_values.GetStrainVector()[d];
            vol_stress += cons_law_values.GetStressVector()[d];
        }
        vol_stress /= dim;
        const double bulk_modulus = (std::abs(vol_strain) > 1.0e-15) ? vol_stress * vol_strain / std::pow(vol_strain, 2) : CalculateApproximatedBulkModulus(rCurrentProcessInfo, i_gauss, kinematic_variables.N);

        // Calculate stabilization constants
        const double h = ComputeElementSize(kinematic_variables.DN_DX);
        const double shear_modulus = cons_law_values.GetConstitutiveMatrix()(dim,dim);
        const double tau_1 = 2.0 * std::pow(h, 2) / (2.0 * shear_modulus);
        const double tau_2 = 0.15;

        // Add the LHS contributions
        const double aux_1 = w_gauss * bulk_modulus * tau_2;
        const double aux_2 = w_gauss * tau_1 * bulk_modulus;
        const double aux_3 = w_gauss * tau_2;
        for (unsigned int i = 0; i < n_nodes; ++i) {
            const unsigned int i_mass_row = i * block_size + dim;
            for (unsigned int j = 0; j < n_nodes; ++j) {
                const unsigned int j_mass_col = j * block_size + dim;
                // Add mass conservation volumetric strain contribution
                rLeftHandSideMatrix(i_mass_row, j_mass_col) -= w_gauss * kinematic_variables.N[i] * kinematic_variables.N[j];
                // Add the volumetric strain stabilization term
                rLeftHandSideMatrix(i_mass_row, j_mass_col) -= aux_3 * kinematic_variables.N(i) * kinematic_variables.N(j);
                for (unsigned int d = 0; d < dim; ++d) {
                    const unsigned int i_mom_row = i * block_size + d;
                    const unsigned int j_mom_col = j * block_size + d;
                    // Add momentum volume stress contribution
                    rLeftHandSideMatrix(i_mom_row, j_mass_col) += w_gauss * bulk_modulus * kinematic_variables.DN_DX(i, d) * kinematic_variables.N[j];
                    // Add mass conservation divergence contribution
                    rLeftHandSideMatrix(i_mass_row, j_mom_col) += w_gauss * kinematic_variables.N[i] * kinematic_variables.DN_DX(j,d);
                    for (unsigned int d2 = 0; d2 < dim; ++d2) {
                        const unsigned int j_mom_col_d2 = j * block_size + d2;
                        // Add the volumetric strain momentum stabilization term - term 1
                        rLeftHandSideMatrix(i_mom_row, j_mom_col_d2) -= aux_1 * kinematic_variables.DN_DX(i, d) * kinematic_variables.DN_DX(j, d2);
                        // Add momentum deviatoric stress contribution
                        for (unsigned int l = 0; l < strain_size; ++l) {
                            for (unsigned int m = 0; m < strain_size; ++m) {
                                for (unsigned int n = 0; n < strain_size; ++n) {
                                    rLeftHandSideMatrix(i_mom_row, j_mom_col_d2) += w_gauss * kinematic_variables.B(l, i * dim + d) * cons_law_values.GetConstitutiveMatrix()(l, m) * kinematic_variables.DevStrainOp(m, n) * kinematic_variables.B(n, j * dim + d2);
                                }
                            }
                        }
                    }
                    // Add the volumetric strain momentum stabilization term - term 2
                    rLeftHandSideMatrix(i_mom_row, j_mass_col) += aux_1 * kinematic_variables.DN_DX(i,d) * kinematic_variables.N(j);
                    // Add the volumetric strain mass stabilization term - term 2
                    rLeftHandSideMatrix(i_mass_row, j_mass_col) -= aux_2 * kinematic_variables.DN_DX(i, d) * kinematic_variables.DN_DX(j, d);
                    // Add the divergence mass stabilization term
                    rLeftHandSideMatrix(i_mass_row, j_mom_col) += aux_3 * kinematic_variables.N(i) * kinematic_variables.DN_DX(j,d);
                }
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    const auto &r_geometry = GetGeometry();
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
    KinematicVariables kinematic_variables(strain_size, dim, n_nodes);
    for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
        const auto &r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        for (unsigned int d = 0; d < dim; ++d) {
            kinematic_variables.Displacements(i_node * dim + d) = r_disp[d];
        }
        kinematic_variables.VolumetricNodalStrains[i_node] = r_geometry[i_node].FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
    }

    // Create the constitutive variables and values containers
    ConstitutiveVariables constitutive_variables(strain_size);
    ConstitutiveLaw::Parameters cons_law_values(r_geometry, GetProperties(), rCurrentProcessInfo);
    auto &r_cons_law_options = cons_law_values.GetOptions();
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // Calculate the RHS contributions
    rRightHandSideVector.clear();

    const SizeType n_gauss = r_geometry.IntegrationPointsNumber(GetIntegrationMethod());
    const auto &r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());
    for (unsigned int i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
        // Calculate kinematics
        CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());
        const double w_gauss = kinematic_variables.detJ0 * r_integration_points[i_gauss].Weight();

        // Calculate the constitutive response
        CalculateConstitutiveVariables(kinematic_variables, constitutive_variables, cons_law_values, i_gauss, r_geometry.IntegrationPoints(this->GetIntegrationMethod()), ConstitutiveLaw::StressMeasure_Cauchy);

        // Calculate stabilization constants
        double vol_strain = 0.0;
        double vol_stress = 0.0;
        for (unsigned int d = 0; d < dim; ++d) {
            vol_strain += cons_law_values.GetStrainVector()[d];
            vol_stress += cons_law_values.GetStressVector()[d];
        }
        vol_stress /= dim;
        const double bulk_modulus = (std::abs(vol_strain) > 1.0e-15) ? vol_stress * vol_strain / std::pow(vol_strain, 2) : CalculateApproximatedBulkModulus(rCurrentProcessInfo, i_gauss, kinematic_variables.N);

        const double h = ComputeElementSize(kinematic_variables.DN_DX);
        const double shear_modulus = cons_law_values.GetConstitutiveMatrix()(dim,dim);
        const double tau_1 = 2.0 * std::pow(h, 2) / (2.0 * shear_modulus);
        const double tau_2 = 0.15;

        // Add the RHS contributions
        const double aux_1 = w_gauss * bulk_modulus * tau_2;
        const double aux_2 = w_gauss * tau_1 * bulk_modulus;
        const double aux_3 = w_gauss * tau_2;
        const auto body_force = GetBodyForce(r_geometry.IntegrationPoints(GetIntegrationMethod()), i_gauss);
        for (unsigned int i = 0; i < n_nodes; ++i) {
            for (unsigned int j = 0; j < n_nodes; ++j) {
                double &r_rhs_mass_row = rRightHandSideVector[i * block_size + dim];
                // Add mass conservation volumetric strain contribution
                r_rhs_mass_row += w_gauss * kinematic_variables.N[i] * kinematic_variables.N[j] * kinematic_variables.VolumetricNodalStrains[j];
                // Add the volumetric strain mass stabilization term - term 2
                r_rhs_mass_row += aux_3 * kinematic_variables.N(i) * kinematic_variables.N(j) * kinematic_variables.VolumetricNodalStrains[j];
                const double aux_mass = kinematic_variables.N[i] * kinematic_variables.N[j];
                for (unsigned int d = 0; d < dim; ++d) {
                    double &r_rhs_mom_row = rRightHandSideVector[i * block_size + d];
                    // Add momentum body force contribution
                    r_rhs_mom_row += w_gauss * aux_mass * body_force[d];
                    // Add momentum stress contribution
                    // Note that this includes both the deviatoric and volumetric stress contributions
                    r_rhs_mom_row -= w_gauss * kinematic_variables.B(j, i * dim + d) * cons_law_values.GetStressVector()[j];
                    // Add mass conservation divergence contribution
                    r_rhs_mass_row -= w_gauss * kinematic_variables.N[i] * kinematic_variables.DN_DX(j, d) * kinematic_variables.Displacements[j * dim + d];
                    // Add the volumetric strain momentum stabilization term - term 2
                    r_rhs_mom_row -= aux_1 * kinematic_variables.DN_DX(i, d) * kinematic_variables.N(j) * kinematic_variables.VolumetricNodalStrains[j];
                    // Add the divergence mass stabilization term - term 2
                    r_rhs_mass_row += aux_2 * kinematic_variables.DN_DX(i, d) * kinematic_variables.DN_DX(j, d) * kinematic_variables.VolumetricNodalStrains[j];
                    // Add the volumetric strain mass stabilization term - term 1
                    r_rhs_mass_row -= aux_3 * kinematic_variables.N(i) * kinematic_variables.DN_DX(j, d) * kinematic_variables.Displacements[j * dim + d];
                    // Add the volumetric strain momentum stabilization term - term 1
                    for (unsigned int d2 = 0; d2 < dim; ++d2) {
                        r_rhs_mom_row += aux_1 * kinematic_variables.DN_DX(i, d) * kinematic_variables.DN_DX(j, d2) * kinematic_variables.Displacements[j * dim + d2];
                    }
                }
            }
        }

        // Add the divergence mass stabilization term - term 1
        for (unsigned int i = 0; i < n_nodes; ++i) {
            double &r_rhs_mass_row = rRightHandSideVector[i * block_size + dim];
            for (unsigned int d = 0; d < dim; ++d) {
                r_rhs_mass_row += w_gauss * tau_1 * kinematic_variables.DN_DX(i,d) * body_force[d];
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::InitializeMaterial()
{
    KRATOS_TRY

    const auto &r_properties = GetProperties();
    if (r_properties[CONSTITUTIVE_LAW] != nullptr) {
        const auto &r_geometry = GetGeometry();
        const auto &r_N_values = r_geometry.ShapeFunctionsValues(mThisIntegrationMethod);
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

bool SmallDisplacementMixedVolumetricStrainElement::UseElementProvidedStrain() const
{
    return true;
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::SetConstitutiveVariables(
    KinematicVariables& rThisKinematicVariables,
    ConstitutiveVariables& rThisConstitutiveVariables,
    ConstitutiveLaw::Parameters& rValues,
    const IndexType PointNumber,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints
    )
{
    // Here we essentially set the input parameters
    rValues.SetStrainVector(rThisKinematicVariables.EquivalentStrain); // equivalent total strain
    rValues.SetShapeFunctionsValues(rThisKinematicVariables.N); // shape functions
    rValues.SetDeterminantF(rThisKinematicVariables.detF); // assuming that det(F) is computed somewhere else
    rValues.SetDeformationGradientF(rThisKinematicVariables.F); // assuming that F is computed somewhere else

    // Here we set the space on which the results shall be written
    rValues.SetConstitutiveMatrix(rThisConstitutiveVariables.D); //assuming the determinant is computed somewhere else
    rValues.SetStressVector(rThisConstitutiveVariables.StressVector); //F computed somewhere else
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::CalculateConstitutiveVariables(
    KinematicVariables& rThisKinematicVariables,
    ConstitutiveVariables& rThisConstitutiveVariables,
    ConstitutiveLaw::Parameters& rValues,
    const IndexType PointNumber,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
    const ConstitutiveLaw::StressMeasure ThisStressMeasure)
{
    // Set the constitutive variables
    SetConstitutiveVariables(rThisKinematicVariables, rThisConstitutiveVariables, rValues, PointNumber, IntegrationPoints);

    // Actually do the computations in the ConstitutiveLaw
    mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(rValues, ThisStressMeasure); //here the calculations are actually done
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double, 3> SmallDisplacementMixedVolumetricStrainElement::GetBodyForce(
    const GeometryType::IntegrationPointsArrayType& rIntegrationPoints,
    const IndexType PointNumber) const
{
    array_1d<double, 3> body_force;
    for (IndexType i = 0; i < 3; ++i) {
        body_force[i] = 0.0;
    }

    const auto &r_properties = GetProperties();
    const double density = r_properties.Has(DENSITY) ? r_properties[DENSITY] : 0.0;

    if (r_properties.Has(VOLUME_ACCELERATION)) {
        noalias(body_force) += density * r_properties[VOLUME_ACCELERATION];
    }

    const auto &r_geometry = GetGeometry();
    if(r_geometry[0].SolutionStepsDataHas(VOLUME_ACCELERATION)) {
        Vector N;
        N = r_geometry.ShapeFunctionsValues(N, rIntegrationPoints[PointNumber].Coordinates());
        for (IndexType i_node = 0; i_node < r_geometry.PointsNumber(); ++i_node) {
            noalias(body_force) += N[i_node] * density * r_geometry[i_node].FastGetSolutionStepValue(VOLUME_ACCELERATION);
        }
    }

    return body_force;
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::CalculateGeometryData(
    const GeometryType &rGeometry,
    Vector &rWeightsContainer,
    Matrix &rShapeFunctionsContainer,
    GeometryType::ShapeFunctionsGradientsType &rDNDX0Container) const
{
    // Compute the geometry data
    const SizeType n_nodes = rGeometry.PointsNumber();
    const auto r_integration_method = GetIntegrationMethod();
    const SizeType n_gauss = rGeometry.IntegrationPointsNumber(r_integration_method);
    const auto& r_integration_points = rGeometry.IntegrationPoints(r_integration_method);

    // Check Gauss points containers size
    if (rShapeFunctionsContainer.size1() != n_gauss || rShapeFunctionsContainer.size2() != n_nodes) {
        rShapeFunctionsContainer.resize(n_gauss, n_nodes,false);
    }
    if (rWeightsContainer.size() != n_gauss) {
        rWeightsContainer.resize(n_gauss, false);
    }
    if (rDNDX0Container.size() != n_gauss){
        rDNDX0Container.resize(n_gauss, false);
    }

    // Calculate the shape function values
    rShapeFunctionsContainer = rGeometry.ShapeFunctionsValues(r_integration_method);

    // Calculate the shape function local gradients
    const auto &rDN_De_container = rGeometry.ShapeFunctionsLocalGradients(r_integration_method);

    for (unsigned int i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
        // Calculate InvJ0
        double detJ0;
        Matrix J0, invJ0;
        GeometryUtils::JacobianOnInitialConfiguration(rGeometry, r_integration_points[i_gauss], J0);
        MathUtils<double>::InvertMatrix(J0, invJ0, detJ0);
        // Calculate shape functions gradients
        GeometryUtils::ShapeFunctionsGradients(rDN_De_container[i_gauss], invJ0, rDNDX0Container[i_gauss]);
        // Calculate integration weights (already multiplied by det(J))
        rWeightsContainer[i_gauss] = detJ0 * r_integration_points[i_gauss].Weight();
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::CalculateKinematicVariables(
    KinematicVariables& rThisKinematicVariables,
    const IndexType PointNumber,
    const GeometryType::IntegrationMethod& rIntegrationMethod) const
{
    const auto &r_geometry = GetGeometry();
    const auto &r_integration_points = r_geometry.IntegrationPoints(rIntegrationMethod);

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

    // Compute B
    CalculateB(rThisKinematicVariables.B, rThisKinematicVariables.DN_DX);

    // Calculate the equivalent total strain
    CalculateEquivalentStrain(rThisKinematicVariables);

    // Compute equivalent F
    ComputeEquivalentF(rThisKinematicVariables.F, rThisKinematicVariables.EquivalentStrain);
    rThisKinematicVariables.detF = MathUtils<double>::Det(rThisKinematicVariables.F);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::CalculateB(
    Matrix& rB,
    const Matrix& rDN_DX) const
{
    KRATOS_TRY;

    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    rB.clear();

    if(dimension == 2) {
        for ( SizeType i = 0; i < number_of_nodes; ++i ) {
            rB(0, i*2    ) = rDN_DX(i, 0);
            rB(1, i*2 + 1) = rDN_DX(i, 1);
            rB(2, i*2    ) = rDN_DX(i, 1);
            rB(2, i*2 + 1) = rDN_DX(i, 0);
        }
    } else if(dimension == 3) {
        for ( SizeType i = 0; i < number_of_nodes; ++i ) {
            rB(0, i*3    ) = rDN_DX(i, 0);
            rB(1, i*3 + 1) = rDN_DX(i, 1);
            rB(2, i*3 + 2) = rDN_DX(i, 2);
            rB(3, i*3    ) = rDN_DX(i, 1);
            rB(3, i*3 + 1) = rDN_DX(i, 0);
            rB(4, i*3 + 1) = rDN_DX(i, 2);
            rB(4, i*3 + 2) = rDN_DX(i, 1);
            rB(5, i*3    ) = rDN_DX(i, 2);
            rB(5, i*3 + 2) = rDN_DX(i, 0);
        }
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::CalculateEquivalentStrain(
    const Vector &rN,
    const Matrix &rB,
    const Matrix &rDevStrainOp,
    Vector &rEquivalentStrain) const
{
    const auto &r_geometry = GetGeometry();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();

    // Initialize the equivalent total strain vector
    rEquivalentStrain = ZeroVector(strain_size);

    // Calculate the equivalent total strain
    // Add the deviatoric contribution
    for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
        const auto &r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        for (unsigned int d = 0; d < dim; ++d) {
            const unsigned int aux = i_node * dim + d;
            for (unsigned int i = 0; i < strain_size; ++i) {
                for (unsigned int j = 0; j < strain_size ; ++j) {
                    rEquivalentStrain[i] += rDevStrainOp(i,j) * rB(j,aux) * r_disp[d];
                }
            }
        }
    }

    // Interpolate and add the nodal volumetric strain
    double gauss_vol_strain = 0.0;
    for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
        gauss_vol_strain += rN[i_node] * r_geometry[i_node].FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
    }
    for (unsigned int d = 0; d < dim; ++d) {
        rEquivalentStrain[d] += (1.0/dim) * gauss_vol_strain;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::CalculateEquivalentStrain(KinematicVariables& rThisKinematicVariables) const
{
    const auto &r_geom = GetGeometry();
    const SizeType n_nodes = r_geom.PointsNumber();
    const SizeType dim = r_geom.WorkingSpaceDimension();

    // Add the deviatoric contribution to the equivalent strain
    const Vector total_strain = prod(rThisKinematicVariables.B, rThisKinematicVariables.Displacements);
    rThisKinematicVariables.EquivalentStrain = prod(rThisKinematicVariables.DevStrainOp, total_strain);

    // Interpolate and add the nodal volumetric strain contribution
    double gauss_vol_strain = 0.0;
    for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
        gauss_vol_strain += rThisKinematicVariables.N[i_node] * rThisKinematicVariables.VolumetricNodalStrains[i_node];
    }
    for (unsigned int d = 0; d < dim; ++d) {
        rThisKinematicVariables.EquivalentStrain[d] += (1.0/dim) * gauss_vol_strain;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::ComputeEquivalentF(
    Matrix &rF,
    const Vector &rStrainTensor) const
{
    const SizeType dim = GetGeometry().WorkingSpaceDimension();

    if(dim == 2) {
        rF(0,0) = 1.0+rStrainTensor(0);
        rF(0,1) = 0.5*rStrainTensor(2);
        rF(1,0) = 0.5*rStrainTensor(2);
        rF(1,1) = 1.0+rStrainTensor(1);
    } else {
        rF(0,0) = 1.0+rStrainTensor(0);
        rF(0,1) = 0.5*rStrainTensor(3);
        rF(0,2) = 0.5*rStrainTensor(5);
        rF(1,0) = 0.5*rStrainTensor(3);
        rF(1,1) = 1.0+rStrainTensor(1);
        rF(1,2) = 0.5*rStrainTensor(4);
        rF(2,0) = 0.5*rStrainTensor(5);
        rF(2,1) = 0.5*rStrainTensor(4);
        rF(2,2) = 1.0+rStrainTensor(2);
    }
}

/***********************************************************************************/
/***********************************************************************************/

double SmallDisplacementMixedVolumetricStrainElement::ComputeElementSize(const Matrix &rDN_DX) const
{
    double h = 0.0;
    for (unsigned int i_node = 0; i_node < GetGeometry().PointsNumber(); ++i_node) {
        double h_inv = 0.0;
        for (unsigned int d = 0; d < GetGeometry().WorkingSpaceDimension(); ++d) {
            h_inv += rDN_DX(i_node, d) * rDN_DX(i_node, d);
        }
        h += 1.0 / h_inv;
    }
    h = sqrt(h) / static_cast<double>(GetGeometry().PointsNumber());
    return h;
}

/***********************************************************************************/
/***********************************************************************************/

double SmallDisplacementMixedVolumetricStrainElement::CalculateApproximatedBulkModulus(
    const ProcessInfo& rCurrentProcessInfo,
    const SizeType i_gauss,
    const Vector &rN) const
{
    const auto &r_geom = GetGeometry();
    const SizeType dim = r_geom.WorkingSpaceDimension();

    // Calculate the bulk modulus with a fake volumetric strain field
    SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();
    Vector aux_vol_stress_vect(strain_size);
    Vector aux_vol_strain_vect = ZeroVector(strain_size);
    for (unsigned int d = 0; d < dim; ++d) {
        aux_vol_strain_vect[d] = 1.0;
    }

    // Call the constitutive law to get the material response of the fake volumetric strain field
    Matrix deformation_gradient(dim, dim);
    ConstitutiveLaw::Parameters aux_cons_law_values(r_geom, GetProperties(), rCurrentProcessInfo);
    auto &r_aux_cons_law_options = aux_cons_law_values.GetOptions();
    r_aux_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_aux_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_aux_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
    aux_cons_law_values.SetShapeFunctionsValues(rN);
    aux_cons_law_values.SetStrainVector(aux_vol_strain_vect);
    aux_cons_law_values.SetStressVector(aux_vol_stress_vect);
    ComputeEquivalentF(deformation_gradient, aux_vol_strain_vect);
    aux_cons_law_values.SetDeformationGradientF(deformation_gradient);
    aux_cons_law_values.SetDeterminantF(MathUtils<double>::Det(deformation_gradient));
    mConstitutiveLawVector[i_gauss]->CalculateMaterialResponseCauchy(aux_cons_law_values);

    double aux_vol_strain = 0.0;
    double aux_vol_stress = 0.0;
    const double alpha = r_geom.WorkingSpaceDimension();
    for (unsigned int d = 0; d < dim; ++d) {
        aux_vol_strain += aux_vol_strain_vect[d];
        aux_vol_stress += aux_cons_law_values.GetStressVector()[d];
    }
    aux_vol_stress /= alpha;

    return aux_vol_stress * aux_vol_strain / std::pow(aux_vol_strain, 2);
}

/***********************************************************************************/
/***********************************************************************************/

int  SmallDisplacementMixedVolumetricStrainElement::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    int base_element_check = SmallDisplacementMixedVolumetricStrainElement::BaseType::Check(rCurrentProcessInfo);

    return base_element_check;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SmallDisplacementMixedVolumetricStrainElement::BaseType);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SmallDisplacementMixedVolumetricStrainElement::BaseType);
}

} // Namespace Kratos
