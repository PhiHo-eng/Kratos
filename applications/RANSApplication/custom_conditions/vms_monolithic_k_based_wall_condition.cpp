//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <string>

// External includes

// Project includes
#include "includes/cfd_variables.h"
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/variables.h"

// Application includes
#include "custom_utilities/fluid_calculation_utilities.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "vms_monolithic_k_based_wall_condition.h"

namespace Kratos
{
template <unsigned int TDim, unsigned int TNumNodes>
int VMSMonolithicKBasedWallCondition<TDim, TNumNodes>::Check(
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    int check = BaseType::Check(rCurrentProcessInfo);

    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(TURBULENCE_RANS_C_MU))
        << "TURBULENCE_RANS_C_MU is not found in process info.\n";
    KRATOS_ERROR_IF(rCurrentProcessInfo[VON_KARMAN] <= 0.0)
        << "VON_KARMAN is not found or not properly initialized in process "
           "info.\n";

    const auto& r_parent_element = this->GetValue(NEIGHBOUR_ELEMENTS)[0];
    const auto& r_element_properties = this->GetValue(NEIGHBOUR_ELEMENTS)[0].GetProperties();
    KRATOS_ERROR_IF(r_element_properties[DENSITY] <= 0.0)
        << "DENSITY is not found or not properly initialized in "
           "properties of "
           "parent element [ rCondition.Id() = "
        << this->Id() << ", rParentElement.Id() = " << r_parent_element.Id()
        << ", rElementProperties.Id() = " << r_element_properties.Id()
        << ", DENSITY = " << r_element_properties[DENSITY] << " ].\n";
    KRATOS_ERROR_IF(r_element_properties[DYNAMIC_VISCOSITY] <= 0.0)
        << "DYNAMIC_VISCOSITY is not found or not properly initialized in "
           "properties of "
           "parent element [ rCondition.Id() = "
        << this->Id() << ", rParentElement.Id() = " << r_parent_element.Id()
        << ", rElementProperties.Id() = " << r_element_properties.Id()
        << ", DYNAMIC_VISCOSITY = " << r_element_properties[DYNAMIC_VISCOSITY] << " ].\n";

    if (RansCalculationUtilities::IsWallFunctionActive(*this)) {

        KRATOS_ERROR_IF_NOT(this->Has(NEIGHBOUR_ELEMENTS))
            << "NEIGHBOUR_ELEMENTS are not found in condition [ rCondition.Id() = "
            << this->Id() << " ].\n";

        KRATOS_ERROR_IF(this->GetValue(NEIGHBOUR_ELEMENTS).size() != 1)
            << "NEIGHBOUR_ELEMENTS are not properly initialized in condition. "
            "Condition can only have one parent [ rCondition.Id() = "
            << this->Id() << ", NEIGHBOUR_ELEMENTS.size() = "
            << this->GetValue(NEIGHBOUR_ELEMENTS).size() << ", required.size() = 1 ].\n";

        KRATOS_ERROR_IF(norm_2(this->GetValue(NORMAL)) == 0.0)
            << "NORMAL is properly initialized in condition. [ rCondition.Id() = "
            << this->Id() << " ].\n";

        const auto& r_condition_properties = this->GetProperties();
        KRATOS_ERROR_IF(r_condition_properties[WALL_SMOOTHNESS_BETA] <= 0.0)
            << "WALL_SMOOTHNESS_BETA is not found or not properly initialized in "
            "properties of "
            "condition [ rCondition.Id() = "
            << this->Id() << ", rConditionProperties.Id() = " << r_condition_properties.Id()
            << ", WALL_SMOOTHNESS_BETA = " << r_condition_properties[WALL_SMOOTHNESS_BETA] << " ].\n";
        KRATOS_ERROR_IF(r_condition_properties[RANS_LINEAR_LOG_LAW_Y_PLUS_LIMIT] <= 0.0)
            << "RANS_LINEAR_LOG_LAW_Y_PLUS_LIMIT is not found or not properly initialized in "
            "properties of "
            "condition [ rCondition.Id() = "
            << this->Id() << ", rConditionProperties.Id() = " << r_condition_properties.Id()
            << ", RANS_LINEAR_LOG_LAW_Y_PLUS_LIMIT = " << r_condition_properties[RANS_LINEAR_LOG_LAW_Y_PLUS_LIMIT] << " ].\n";
    }

    const auto& r_geometry = this->GetGeometry();

    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
        const auto& r_node = r_geometry[i_node];

        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(MESH_VELOCITY, r_node);
    }

    return check;

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void VMSMonolithicKBasedWallCondition<TDim, TNumNodes>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    if (RansCalculationUtilities::IsWallFunctionActive(*this)) {
        const array_1d<double, 3>& r_normal = this->GetValue(NORMAL);

        mWallHeight = RansCalculationUtilities::CalculateWallHeight(*this, r_normal);

        this->SetValue(GAUSS_RANS_Y_PLUS, Vector(this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod())));

        this->SetValue(DISTANCE, mWallHeight);

        KRATOS_ERROR_IF(mWallHeight == 0.0) << this->Info() << " has zero wall height.\n";
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void VMSMonolithicKBasedWallCondition<TDim, TNumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    const SizeType BlockSize = TDim + 1;
    const SizeType LocalSize = BlockSize * TNumNodes;

    if (rLeftHandSideMatrix.size1() != LocalSize)
        rLeftHandSideMatrix.resize(LocalSize,LocalSize);

    if (rRightHandSideVector.size() != LocalSize)
        rRightHandSideVector.resize(LocalSize);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize,LocalSize);
    noalias(rRightHandSideVector) = ZeroVector(LocalSize);
}

template <unsigned int TDim, unsigned int TNumNodes>
void VMSMonolithicKBasedWallCondition<TDim, TNumNodes>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    const SizeType BlockSize = TDim + 1;
    const SizeType LocalSize = BlockSize * TNumNodes;

    if (rLeftHandSideMatrix.size1() != LocalSize)
        rLeftHandSideMatrix.resize(LocalSize,LocalSize);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize,LocalSize);
}

template <unsigned int TDim, unsigned int TNumNodes>
void VMSMonolithicKBasedWallCondition<TDim, TNumNodes>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    const SizeType BlockSize = TDim + 1;
    const SizeType LocalSize = BlockSize * TNumNodes;

    if (rRightHandSideVector.size() != LocalSize)
        rRightHandSideVector.resize(LocalSize);

    noalias(rRightHandSideVector) = ZeroVector(LocalSize);
}

template <unsigned int TDim, unsigned int TNumNodes>
void VMSMonolithicKBasedWallCondition<TDim, TNumNodes>::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    VectorType RHS;
    this->CalculateLocalVelocityContribution(rDampingMatrix,RHS,rCurrentProcessInfo);
}

template<unsigned int TDim, unsigned int TNumNodes>
void VMSMonolithicKBasedWallCondition<TDim, TNumNodes>::CalculateLocalVelocityContribution(
    MatrixType &rDampMatrix,
    VectorType &rRightHandSideVector,
    const ProcessInfo &rCurrentProcessInfo)
{
    // Initialize local contributions
    const SizeType LocalSize = (TDim + 1) * TNumNodes;

    if (rDampMatrix.size1() != LocalSize)
        rDampMatrix.resize(LocalSize,LocalSize);
    if (rRightHandSideVector.size() != LocalSize)
        rRightHandSideVector.resize(LocalSize);

    noalias(rDampMatrix) = ZeroMatrix(LocalSize,LocalSize);
    noalias(rRightHandSideVector) = ZeroVector(LocalSize);

    this->ApplyWallLaw(rDampMatrix,rRightHandSideVector,rCurrentProcessInfo);
}

template <>
void VMSMonolithicKBasedWallCondition<2, 2>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const unsigned int NumNodes = 2;
    const unsigned int LocalSize = 6;
    unsigned int LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_X).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Y).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(PRESSURE).EquationId();
    }
}

template <>
void VMSMonolithicKBasedWallCondition<3, 3>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const SizeType NumNodes = 3;
    const SizeType LocalSize = 12;
    unsigned int LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_X).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Y).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Z).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(PRESSURE).EquationId();
    }
}

template <>
void VMSMonolithicKBasedWallCondition<2, 2>::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const SizeType NumNodes = 2;
    const SizeType LocalSize = 6;

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    unsigned int LocalIndex = 0;

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_X);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Y);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(PRESSURE);
    }
}

template <>
void VMSMonolithicKBasedWallCondition<3, 3>::GetDofList(DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const SizeType NumNodes = 3;
    const SizeType LocalSize = 12;

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    unsigned int LocalIndex = 0;

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_X);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Y);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Z);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(PRESSURE);
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void VMSMonolithicKBasedWallCondition<TDim, TNumNodes>::GetFirstDerivativesVector(
    Vector& Values,
    int Step) const
{
    const SizeType LocalSize = (TDim + 1) * TNumNodes;
    unsigned int LocalIndex = 0;

    if (Values.size() != LocalSize)
        Values.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
    {
        const array_1d<double,3>& rVelocity = this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY, Step);
        for (unsigned int d = 0; d < TDim; ++d)
            Values[LocalIndex++] = rVelocity[d];
        Values[LocalIndex++] = this->GetGeometry()[iNode].FastGetSolutionStepValue(PRESSURE, Step);
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void VMSMonolithicKBasedWallCondition<TDim, TNumNodes>::GetSecondDerivativesVector(
    Vector& Values,
    int Step) const
{
    const SizeType LocalSize = (TDim + 1) * TNumNodes;
    unsigned int LocalIndex = 0;

    if (Values.size() != LocalSize)
        Values.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
    {
        const array_1d<double,3>& rVelocity = this->GetGeometry()[iNode].FastGetSolutionStepValue(ACCELERATION, Step);
        for (unsigned int d = 0; d < TDim; ++d)
            Values[LocalIndex++] = rVelocity[d];
        Values[LocalIndex++] = 0.0; // No value on pressure positions
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
std::string VMSMonolithicKBasedWallCondition<TDim, TNumNodes>::Info() const
{
    std::stringstream buffer;
    buffer << "VMSMonolithicKBasedWallCondition" << TDim << "D";
    return buffer.str();
}

template <unsigned int TDim, unsigned int TNumNodes>
void VMSMonolithicKBasedWallCondition<TDim, TNumNodes>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "VMSMonolithicKBasedWallCondition";
}

template <unsigned int TDim, unsigned int TNumNodes>
void VMSMonolithicKBasedWallCondition<TDim, TNumNodes>::PrintData(std::ostream& rOStream) const
{
}

template <unsigned int TDim, unsigned int TNumNodes>
void VMSMonolithicKBasedWallCondition<TDim, TNumNodes>::ApplyWallLaw(
    MatrixType& rLocalMatrix,
    VectorType& rLocalVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    using namespace RansCalculationUtilities;

    if (IsWallFunctionActive(*this)) {
        const auto& r_geometry = this->GetGeometry();
        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        CalculateConditionGeometryData(r_geometry, this->GetIntegrationMethod(),
                                       gauss_weights, shape_functions);
        const IndexType num_gauss_points = gauss_weights.size();

        const size_t block_size = TDim + 1;

        const double eps = std::numeric_limits<double>::epsilon();
        const double c_mu_25 = std::pow(rCurrentProcessInfo[TURBULENCE_RANS_C_MU], 0.25);
        const double kappa = rCurrentProcessInfo[VON_KARMAN];

        // get parent element
        auto& r_parent_element = this->GetValue(NEIGHBOUR_ELEMENTS)[0];

        // get fluid properties from parent element
        const auto& r_elem_properties = r_parent_element.GetProperties();
        const double rho = r_elem_properties[DENSITY];
        const double nu = r_elem_properties[DYNAMIC_VISCOSITY] / rho;

        // get surface properties from condition
        const PropertiesType& r_cond_properties = this->GetProperties();
        const double beta = r_cond_properties.GetValue(WALL_SMOOTHNESS_BETA);
        const double y_plus_limit = r_cond_properties.GetValue(RANS_LINEAR_LOG_LAW_Y_PLUS_LIMIT);
        const double inv_kappa = 1.0 / kappa;

        double tke;
        array_1d<double, 3> wall_velocity, fluid_velocity, mesh_velocity;

        auto& gauss_y_plus = this->GetValue(GAUSS_RANS_Y_PLUS);

        for (size_t g = 0; g < num_gauss_points; ++g) {
            const Vector& gauss_shape_functions = row(shape_functions, g);
            double& y_plus = gauss_y_plus[g];

            FluidCalculationUtilities::EvaluateInPoint(
                r_geometry, gauss_shape_functions,
                std::tie(tke, TURBULENT_KINETIC_ENERGY),
                std::tie(fluid_velocity, VELOCITY),
                std::tie(mesh_velocity, MESH_VELOCITY));

            noalias(wall_velocity) = fluid_velocity - mesh_velocity;

            const double wall_velocity_magnitude = norm_2(wall_velocity);

            double u_tau{0.0};
            CalculateYPlusAndUtau(y_plus, u_tau, wall_velocity_magnitude,
                                  mWallHeight, nu, kappa, beta);
            y_plus = std::max(y_plus, y_plus_limit);

            if (wall_velocity_magnitude > eps) {
                const double tke_u_tau = c_mu_25 * std::sqrt(std::max(tke, 0.0));
                const double u_u_tau = wall_velocity_magnitude / (inv_kappa * std::log(y_plus) + beta);
                const double u_tau = std::max(tke_u_tau, u_u_tau);

                const double value = rho * std::pow(u_tau, 2) *
                                     gauss_weights[g] / wall_velocity_magnitude;

                for (IndexType a = 0; a < r_geometry.PointsNumber(); ++a) {
                    for (IndexType dim = 0; dim < TDim; ++dim) {
                        for (IndexType b = 0; b < r_geometry.PointsNumber(); ++b) {
                            rLocalMatrix(a * block_size + dim, b * block_size + dim) +=
                                gauss_shape_functions[a] * gauss_shape_functions[b] * value;
                        }
                        rLocalVector[a * block_size + dim] -=
                            gauss_shape_functions[a] * value * wall_velocity[dim];
                    }
                }
            }
        }
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
GeometryData::IntegrationMethod VMSMonolithicKBasedWallCondition<TDim, TNumNodes>::GetIntegrationMethod()
{
    return GeometryData::IntegrationMethod::GI_GAUSS_2;
}

template <unsigned int TDim, unsigned int TNumNodes>
void VMSMonolithicKBasedWallCondition<TDim, TNumNodes>::save(Serializer& rSerializer) const
{
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
}

template <unsigned int TDim, unsigned int TNumNodes>
void VMSMonolithicKBasedWallCondition<TDim, TNumNodes>::load(Serializer& rSerializer)
{
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
}

// template instantiations

template class VMSMonolithicKBasedWallCondition<2, 2>;
template class VMSMonolithicKBasedWallCondition<3, 3>;

} // namespace Kratos.
