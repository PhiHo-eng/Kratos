//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Eduard Gómez
//

// System includes


// External includes


// Project includes
#include "includes/checks.h"
#include "includes/key_hash.h"
#include "includes/ublas_interface.h"
#include "utilities/atomic_utilities.h"
#include "utilities/element_size_calculator.h"

// Application includes
#include "custom_elements/compressible_navier_stokes_explicit.h"


namespace Kratos {

template <>
void CompressibleNavierStokesExplicit<2,4>::EquationIdVector(
    EquationIdVectorType &rResult,
    const ProcessInfo &rCurrentProcessInfo) const
{
    KRATOS_TRY

    if (rResult.size() != DofSize) {
        rResult.resize(DofSize);
    }

    unsigned int local_index = 0;
    const auto& r_geometry = GetGeometry();
    const unsigned int den_pos = r_geometry[0].GetDofPosition(DENSITY);
    const unsigned int mom_pos = r_geometry[0].GetDofPosition(MOMENTUM);
    const unsigned int enr_pos = r_geometry[0].GetDofPosition(TOTAL_ENERGY);
    for (unsigned int i_node = 0; i_node < NumNodes; ++i_node) {
        rResult[local_index++] = r_geometry[i_node].GetDof(DENSITY, den_pos).EquationId();
        rResult[local_index++] = r_geometry[i_node].GetDof(MOMENTUM_X, mom_pos).EquationId();
        rResult[local_index++] = r_geometry[i_node].GetDof(MOMENTUM_Y, mom_pos + 1).EquationId();
        rResult[local_index++] = r_geometry[i_node].GetDof(TOTAL_ENERGY, enr_pos).EquationId();
    }

    KRATOS_CATCH("");
}

template <>
void CompressibleNavierStokesExplicit<2,4>::GetDofList(
    DofsVectorType &ElementalDofList,
    const ProcessInfo &rCurrentProcessInfo) const
{
    KRATOS_TRY

    if (ElementalDofList.size() != DofSize) {
        ElementalDofList.resize(DofSize);
    }

    unsigned int local_index = 0;
    const auto& r_geometry = GetGeometry();
    const unsigned int den_pos = r_geometry[0].GetDofPosition(DENSITY);
    const unsigned int mom_pos = r_geometry[0].GetDofPosition(MOMENTUM);
    const unsigned int enr_pos = r_geometry[0].GetDofPosition(TOTAL_ENERGY);
    for (unsigned int i_node = 0; i_node < NumNodes; ++i_node) {
        ElementalDofList[local_index++] = r_geometry[i_node].pGetDof(DENSITY, den_pos);
        ElementalDofList[local_index++] = r_geometry[i_node].pGetDof(MOMENTUM_X, mom_pos);
        ElementalDofList[local_index++] = r_geometry[i_node].pGetDof(MOMENTUM_Y, mom_pos + 1);
        ElementalDofList[local_index++] = r_geometry[i_node].pGetDof(TOTAL_ENERGY, enr_pos);
    }

    KRATOS_CATCH("");
}


template <>
array_1d<double,3> CompressibleNavierStokesExplicit<2,4>::CalculateMidPointVelocityRotational() const
{
    // Get geometry data
    const auto& r_geom = GetGeometry();
    Geometry<Node<3>>::ShapeFunctionsGradientsType dNdX_container;
    r_geom.ShapeFunctionsIntegrationPointsGradients(dNdX_container, GeometryData::IntegrationMethod::GI_GAUSS_1);
    const auto& r_dNdX = dNdX_container[0];

    // Calculate midpoint magnitudes
    double midpoint_rho = 0.0;
    double midpoint_dmy_dx = 0.0;
    double midpoint_dmx_dy = 0.0;
    double midpoint_rho_dx = 0.0;
    double midpoint_rho_dy = 0.0;
    array_1d<double,3> midpoint_mom = ZeroVector(3);
    for (unsigned int i_node = 0; i_node < NumNodes; ++i_node) {
        auto& r_node = r_geom[i_node];
        const auto node_dNdX = row(r_dNdX, i_node);
        const auto& r_mom = r_node.FastGetSolutionStepValue(MOMENTUM);
        const double& r_rho = r_node.FastGetSolutionStepValue(DENSITY);
        midpoint_rho += r_rho;
        midpoint_mom += r_mom;
        midpoint_dmy_dx += r_mom[1] * node_dNdX[0];
        midpoint_dmx_dy += r_mom[0] * node_dNdX[1];
        midpoint_rho_dx += r_rho * node_dNdX[0];
        midpoint_rho_dy += r_rho * node_dNdX[1];
    }
    midpoint_rho /= NumNodes;
    midpoint_mom /= NumNodes;

    // Calculate velocity rotational
    // Note that the formulation is written in conservative variables. Hence we do rot(mom/rho).
    const double dvy_dx = (midpoint_dmy_dx * midpoint_rho - midpoint_mom[1] * midpoint_rho_dx) / std::pow(midpoint_rho, 2);
    const double dvx_dy = (midpoint_dmx_dy * midpoint_rho - midpoint_mom[0] * midpoint_rho_dy) / std::pow(midpoint_rho, 2);
    array_1d<double,3> midpoint_rot_v;
    midpoint_rot_v[0] = 0.0;
    midpoint_rot_v[1] = 0.0;
    midpoint_rot_v[2] = dvy_dx - dvx_dy;
    return midpoint_rot_v;
}

template <>
BoundedMatrix<double, 3, 3> CompressibleNavierStokesExplicit<2, 4>::CalculateMidPointVelocityGradient() const
{
    KRATOS_TRY

    // Get geometry data
    const auto& r_geom = GetGeometry();
    Geometry<Node<3>>::ShapeFunctionsGradientsType dNdX_container;
    r_geom.ShapeFunctionsIntegrationPointsGradients(dNdX_container, GeometryData::IntegrationMethod::GI_GAUSS_1);
    const auto& r_dNdX = dNdX_container[0];

    // Calculate midpoint magnitudes
    double midpoint_rho = 0.0;
    double midpoint_dmx_dx = 0.0;
    double midpoint_dmx_dy = 0.0;
    double midpoint_dmy_dx = 0.0;
    double midpoint_dmy_dy = 0.0;
    double midpoint_rho_dx = 0.0;
    double midpoint_rho_dy = 0.0;
    array_1d<double,3> midpoint_mom = ZeroVector(3);
    for (unsigned int i_node = 0; i_node < NumNodes; ++i_node) {
        auto& r_node = r_geom[i_node];
        const auto node_dNdX = row(r_dNdX, i_node);
        const auto& r_mom = r_node.FastGetSolutionStepValue(MOMENTUM);
        const double& r_rho = r_node.FastGetSolutionStepValue(DENSITY);
        midpoint_rho += r_rho;
        midpoint_mom += r_mom;
        midpoint_dmx_dx += r_mom[0] * node_dNdX[0];
        midpoint_dmx_dy += r_mom[0] * node_dNdX[1];
        midpoint_dmy_dx += r_mom[1] * node_dNdX[0];
        midpoint_dmy_dy += r_mom[1] * node_dNdX[1];
        midpoint_rho_dx += r_rho * node_dNdX[0];
        midpoint_rho_dy += r_rho * node_dNdX[1];
    }
    midpoint_rho /= NumNodes;
    midpoint_mom /= NumNodes;

    // Calculate velocity gradient
    // Note that the formulation is written in conservative variables. Hence we do grad(mom/rho).
    BoundedMatrix<double, 3, 3> midpoint_grad_v = ZeroMatrix(3, 3);
    midpoint_grad_v(0,0) = (midpoint_dmx_dx * midpoint_rho - midpoint_mom[0] * midpoint_rho_dx);
    midpoint_grad_v(0,1) = (midpoint_dmx_dy * midpoint_rho - midpoint_mom[0] * midpoint_rho_dy);
    midpoint_grad_v(1,0) = (midpoint_dmy_dx * midpoint_rho - midpoint_mom[1] * midpoint_rho_dx);
    midpoint_grad_v(1,1) = (midpoint_dmy_dy * midpoint_rho - midpoint_mom[1] * midpoint_rho_dy);
    midpoint_grad_v /= std::pow(midpoint_rho, 2);

    return midpoint_grad_v;

    KRATOS_CATCH("")
}

template <>
void CompressibleNavierStokesExplicit<2,4>::CalculateMomentumProjection(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Calculate shock capturing values
    BoundedVector<double, Dim*NumNodes> mom_proj = ZeroVector(Dim*NumNodes);

    Vector N;
    Matrix DN_DX_iso;
    Matrix DN_DX;
    Matrix Jinv;

    auto& r_geometry = GetGeometry();
    const auto& gauss_points = r_geometry.IntegrationPoints(GeometryData::IntegrationMethod::GI_GAUSS_2);
    for(const auto& gauss_point: gauss_points)
    {
        r_geometry.ShapeFunctionsValues(N, gauss_point.Coordinates());
        r_geometry.InverseOfJacobian(Jinv, gauss_point.Coordinates());
        r_geometry.ShapeFunctionsLocalGradients(DN_DX_iso, gauss_point.Coordinates());
        GeometryUtils::ShapeFunctionsGradients(DN_DX_iso, Jinv, DN_DX);

const double cmom_proj0 =             data.gamma - 1;
const double cmom_proj1 =             N(0)*data.U(0,0) + N(1)*data.U(1,0) + N(2)*data.U(2,0) + N(3)*data.U(3,0);
const double cmom_proj2 =             DN_DX(0,1)*data.U(0,1) + DN_DX(1,1)*data.U(1,1) + DN_DX(2,1)*data.U(2,1) + DN_DX(3,1)*data.U(3,1);
const double cmom_proj3 =             1.0/cmom_proj1;
const double cmom_proj4 =             N(0)*data.U(0,2) + N(1)*data.U(1,2) + N(2)*data.U(2,2) + N(3)*data.U(3,2);
const double cmom_proj5 =             cmom_proj3*cmom_proj4;
const double cmom_proj6 =             DN_DX(0,1)*data.U(0,2) + DN_DX(1,1)*data.U(1,2) + DN_DX(2,1)*data.U(2,2) + DN_DX(3,1)*data.U(3,2);
const double cmom_proj7 =             N(0)*data.U(0,1) + N(1)*data.U(1,1) + N(2)*data.U(2,1) + N(3)*data.U(3,1);
const double cmom_proj8 =             cmom_proj3*cmom_proj7;
const double cmom_proj9 =             DN_DX(0,0)*data.U(0,2) + DN_DX(1,0)*data.U(1,2) + DN_DX(2,0)*data.U(2,2) + DN_DX(3,0)*data.U(3,2);
const double cmom_proj10 =             1.0*cmom_proj0;
const double cmom_proj11 =             DN_DX(0,0)*data.U(0,1) + DN_DX(1,0)*data.U(1,1) + DN_DX(2,0)*data.U(2,1) + DN_DX(3,0)*data.U(3,1);
const double cmom_proj12 =             1.0*data.gamma - 3.0;
const double cmom_proj13 =             pow(cmom_proj1, -2);
const double cmom_proj14 =             cmom_proj13*(DN_DX(0,1)*data.U(0,0) + DN_DX(1,1)*data.U(1,0) + DN_DX(2,1)*data.U(2,0) + DN_DX(3,1)*data.U(3,0));
const double cmom_proj15 =             cmom_proj4*cmom_proj7;
const double cmom_proj16 =             pow(cmom_proj7, 2);
const double cmom_proj17 =             pow(cmom_proj4, 2);
const double cmom_proj18 =             0.5*cmom_proj0*(cmom_proj16 + cmom_proj17);
const double cmom_proj19 =             cmom_proj13*(DN_DX(0,0)*data.U(0,0) + DN_DX(1,0)*data.U(1,0) + DN_DX(2,0)*data.U(2,0) + DN_DX(3,0)*data.U(3,0));
const double cmom_proj20 =             N(0)*data.dUdt(0,1) + N(1)*data.dUdt(1,1) + N(2)*data.dUdt(2,1) + N(3)*data.dUdt(3,1) + cmom_proj0*(DN_DX(0,0)*data.U(0,3) + DN_DX(1,0)*data.U(1,3) + DN_DX(2,0)*data.U(2,3) + DN_DX(3,0)*data.U(3,3)) - cmom_proj1*(N(0)*data.f_ext(0,0) + N(1)*data.f_ext(1,0) + N(2)*data.f_ext(2,0) + N(3)*data.f_ext(3,0)) - cmom_proj10*cmom_proj5*cmom_proj9 - cmom_proj11*cmom_proj12*cmom_proj8 - cmom_proj14*cmom_proj15 + cmom_proj19*(-cmom_proj16 + cmom_proj18) + cmom_proj2*cmom_proj5 + cmom_proj6*cmom_proj8;
const double cmom_proj21 =             N(0)*data.dUdt(0,2) + N(1)*data.dUdt(1,2) + N(2)*data.dUdt(2,2) + N(3)*data.dUdt(3,2) + cmom_proj0*(DN_DX(0,1)*data.U(0,3) + DN_DX(1,1)*data.U(1,3) + DN_DX(2,1)*data.U(2,3) + DN_DX(3,1)*data.U(3,3)) - cmom_proj1*(N(0)*data.f_ext(0,1) + N(1)*data.f_ext(1,1) + N(2)*data.f_ext(2,1) + N(3)*data.f_ext(3,1)) - cmom_proj10*cmom_proj2*cmom_proj8 + cmom_proj11*cmom_proj5 - cmom_proj12*cmom_proj5*cmom_proj6 + cmom_proj14*(-cmom_proj17 + cmom_proj18) - cmom_proj15*cmom_proj19 + cmom_proj8*cmom_proj9;
            mom_proj[0] += -N(0)*cmom_proj20;
            mom_proj[1] += -N(0)*cmom_proj21;
            mom_proj[2] += -N(1)*cmom_proj20;
            mom_proj[3] += -N(1)*cmom_proj21;
            mom_proj[4] += -N(2)*cmom_proj20;
            mom_proj[5] += -N(2)*cmom_proj21;
            mom_proj[6] += -N(3)*cmom_proj20;
            mom_proj[7] += -N(3)*cmom_proj21;

    }

    // Here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/NumNodes
    mom_proj *= data.volume / NumNodes;

    // Assembly the projection contributions
    for (IndexType i_node = 0; i_node < NumNodes; ++i_node) {
        const IndexType aux = i_node * Dim;
        auto& r_mom_proj = r_geometry[i_node].GetValue(MOMENTUM_PROJECTION);
        for (IndexType d = 0; d < Dim; ++d) {
            AtomicAdd(r_mom_proj[d], mom_proj[aux + d]);
        }
    }

    KRATOS_CATCH("")
}


template <>
void CompressibleNavierStokesExplicit<2,4>::CalculateDensityProjection(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Calculate shock capturing values
    BoundedVector<double, NumNodes> rho_proj = ZeroVector(NumNodes);

    Vector N;
    Matrix DN_DX_iso;
    Matrix DN_DX;
    Matrix Jinv;

    auto& r_geometry = GetGeometry();
    const auto& gauss_points = r_geometry.IntegrationPoints(GeometryData::IntegrationMethod::GI_GAUSS_2);
    for(const auto& gauss_point: gauss_points)
    {
        r_geometry.ShapeFunctionsValues(N, gauss_point.Coordinates());
        r_geometry.InverseOfJacobian(Jinv, gauss_point.Coordinates());
        r_geometry.ShapeFunctionsLocalGradients(DN_DX_iso, gauss_point.Coordinates());
        GeometryUtils::ShapeFunctionsGradients(DN_DX_iso, Jinv, DN_DX);

const double crho_proj0 =             DN_DX(0,0)*data.U(0,1) + DN_DX(0,1)*data.U(0,2) + DN_DX(1,0)*data.U(1,1) + DN_DX(1,1)*data.U(1,2) + DN_DX(2,0)*data.U(2,1) + DN_DX(2,1)*data.U(2,2) + DN_DX(3,0)*data.U(3,1) + DN_DX(3,1)*data.U(3,2) + N(0)*data.dUdt(0,0) - N(0)*data.m_ext(0) + N(1)*data.dUdt(1,0) - N(1)*data.m_ext(1) + N(2)*data.dUdt(2,0) - N(2)*data.m_ext(2) + N(3)*data.dUdt(3,0) - N(3)*data.m_ext(3);
            rho_proj[0] += -N(0)*crho_proj0;
            rho_proj[1] += -N(1)*crho_proj0;
            rho_proj[2] += -N(2)*crho_proj0;
            rho_proj[3] += -N(3)*crho_proj0;

    }

    // Here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/NumNodes
    rho_proj *= data.volume / NumNodes;

    // Assembly the projection contributions
    for (IndexType i_node = 0; i_node < NumNodes; ++i_node) {
        AtomicAdd(r_geometry[i_node].GetValue(DENSITY_PROJECTION), rho_proj[i_node]);
    }

    KRATOS_CATCH("")
}


template <>
void CompressibleNavierStokesExplicit<2,4>::CalculateTotalEnergyProjection(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Calculate shock capturing values
    BoundedVector<double, NumNodes> tot_ener_proj = ZeroVector(NumNodes);

    Vector N;
    Matrix DN_DX_iso;
    Matrix DN_DX;
    Matrix Jinv;

    auto& r_geometry = GetGeometry();
    const auto& gauss_points = r_geometry.IntegrationPoints(GeometryData::IntegrationMethod::GI_GAUSS_2);
    for(const auto& gauss_point: gauss_points)
    {
        r_geometry.ShapeFunctionsValues(N, gauss_point.Coordinates());
        r_geometry.InverseOfJacobian(Jinv, gauss_point.Coordinates());
        r_geometry.ShapeFunctionsLocalGradients(DN_DX_iso, gauss_point.Coordinates());
        GeometryUtils::ShapeFunctionsGradients(DN_DX_iso, Jinv, DN_DX);

const double ctot_ener_proj0 =             N(0)*data.U(0,0) + N(1)*data.U(1,0) + N(2)*data.U(2,0) + N(3)*data.U(3,0);
const double ctot_ener_proj1 =             N(0)*data.U(0,1) + N(1)*data.U(1,1) + N(2)*data.U(2,1) + N(3)*data.U(3,1);
const double ctot_ener_proj2 =             N(0)*data.U(0,2) + N(1)*data.U(1,2) + N(2)*data.U(2,2) + N(3)*data.U(3,2);
const double ctot_ener_proj3 =             1.0/ctot_ener_proj0;
const double ctot_ener_proj4 =             ctot_ener_proj3*data.gamma;
const double ctot_ener_proj5 =             data.gamma - 1;
const double ctot_ener_proj6 =             pow(ctot_ener_proj0, -2);
const double ctot_ener_proj7 =             1.0*ctot_ener_proj1*ctot_ener_proj2*ctot_ener_proj5*ctot_ener_proj6;
const double ctot_ener_proj8 =             pow(ctot_ener_proj1, 2);
const double ctot_ener_proj9 =             ctot_ener_proj3*ctot_ener_proj5;
const double ctot_ener_proj10 =             1.0*ctot_ener_proj9;
const double ctot_ener_proj11 =             N(0)*data.U(0,3);
const double ctot_ener_proj12 =             N(1)*data.U(1,3);
const double ctot_ener_proj13 =             N(2)*data.U(2,3);
const double ctot_ener_proj14 =             N(3)*data.U(3,3);
const double ctot_ener_proj15 =             pow(ctot_ener_proj2, 2);
const double ctot_ener_proj16 =             -ctot_ener_proj11 - ctot_ener_proj12 - ctot_ener_proj13 - ctot_ener_proj14 - ctot_ener_proj5*(ctot_ener_proj11 + ctot_ener_proj12 + ctot_ener_proj13 + ctot_ener_proj14 - ctot_ener_proj3*(0.5*ctot_ener_proj15 + 0.5*ctot_ener_proj8));
const double ctot_ener_proj17 =             ctot_ener_proj6*(ctot_ener_proj16 + 0.5*ctot_ener_proj9*(ctot_ener_proj15 + ctot_ener_proj8));
const double ctot_ener_proj18 =             N(0)*data.dUdt(0,3) + N(1)*data.dUdt(1,3) + N(2)*data.dUdt(2,3) + N(3)*data.dUdt(3,3) - ctot_ener_proj0*(N(0)*data.r_ext(0) + N(1)*data.r_ext(1) + N(2)*data.r_ext(2) + N(3)*data.r_ext(3)) + ctot_ener_proj1*ctot_ener_proj17*(DN_DX(0,0)*data.U(0,0) + DN_DX(1,0)*data.U(1,0) + DN_DX(2,0)*data.U(2,0) + DN_DX(3,0)*data.U(3,0)) + ctot_ener_proj1*ctot_ener_proj4*(DN_DX(0,0)*data.U(0,3) + DN_DX(1,0)*data.U(1,3) + DN_DX(2,0)*data.U(2,3) + DN_DX(3,0)*data.U(3,3)) - ctot_ener_proj1*(N(0)*data.f_ext(0,0) + N(1)*data.f_ext(1,0) + N(2)*data.f_ext(2,0) + N(3)*data.f_ext(3,0)) + ctot_ener_proj17*ctot_ener_proj2*(DN_DX(0,1)*data.U(0,0) + DN_DX(1,1)*data.U(1,0) + DN_DX(2,1)*data.U(2,0) + DN_DX(3,1)*data.U(3,0)) + ctot_ener_proj2*ctot_ener_proj4*(DN_DX(0,1)*data.U(0,3) + DN_DX(1,1)*data.U(1,3) + DN_DX(2,1)*data.U(2,3) + DN_DX(3,1)*data.U(3,3)) - ctot_ener_proj2*(N(0)*data.f_ext(0,1) + N(1)*data.f_ext(1,1) + N(2)*data.f_ext(2,1) + N(3)*data.f_ext(3,1)) - ctot_ener_proj3*(ctot_ener_proj10*ctot_ener_proj15 + ctot_ener_proj16)*(DN_DX(0,1)*data.U(0,2) + DN_DX(1,1)*data.U(1,2) + DN_DX(2,1)*data.U(2,2) + DN_DX(3,1)*data.U(3,2)) - ctot_ener_proj3*(ctot_ener_proj10*ctot_ener_proj8 + ctot_ener_proj16)*(DN_DX(0,0)*data.U(0,1) + DN_DX(1,0)*data.U(1,1) + DN_DX(2,0)*data.U(2,1) + DN_DX(3,0)*data.U(3,1)) - ctot_ener_proj7*(DN_DX(0,0)*data.U(0,2) + DN_DX(1,0)*data.U(1,2) + DN_DX(2,0)*data.U(2,2) + DN_DX(3,0)*data.U(3,2)) - ctot_ener_proj7*(DN_DX(0,1)*data.U(0,1) + DN_DX(1,1)*data.U(1,1) + DN_DX(2,1)*data.U(2,1) + DN_DX(3,1)*data.U(3,1));
            tot_ener_proj[0] += -N(0)*ctot_ener_proj18;
            tot_ener_proj[1] += -N(1)*ctot_ener_proj18;
            tot_ener_proj[2] += -N(2)*ctot_ener_proj18;
            tot_ener_proj[3] += -N(3)*ctot_ener_proj18;

    }

    // Here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/NumNodes
    tot_ener_proj *= data.volume / NumNodes;

    // Assembly the projection contributions
    for (IndexType i_node = 0; i_node < NumNodes; ++i_node) {
        AtomicAdd(r_geometry[i_node].GetValue(TOTAL_ENERGY_PROJECTION), tot_ener_proj[i_node]);
    }

    KRATOS_CATCH("")
}


template <>
void CompressibleNavierStokesExplicit<2,4>::CalculateRightHandSideInternal(
    BoundedVector<double, DofSize> &rRightHandSideBoundedVector,
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY

    std::fill(rRightHandSideBoundedVector.begin(), rRightHandSideBoundedVector.end(), 0.0);

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Stabilization parameters
    constexpr double stab_c1 = 12.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 1.0;

    const auto& r_geometry = GetGeometry();
    const auto& gauss_points = r_geometry.IntegrationPoints(GeometryData::IntegrationMethod::GI_GAUSS_2);

    Vector N;
    Matrix DN_DX_iso;
    Matrix DN_DX;
    Matrix Jinv;

    if (data.UseOSS)
    {
        for(const auto& gauss_point: gauss_points)
        {
            r_geometry.ShapeFunctionsValues(N, gauss_point.Coordinates());
            r_geometry.InverseOfJacobian(Jinv, gauss_point.Coordinates());
            r_geometry.ShapeFunctionsLocalGradients(DN_DX_iso, gauss_point.Coordinates());
            GeometryUtils::ShapeFunctionsGradients(DN_DX_iso, Jinv, DN_DX);

const double crRightHandSideBoundedVector0 =             N(0)*data.m_ext(0);
const double crRightHandSideBoundedVector1 =             N(1)*data.m_ext(1);
const double crRightHandSideBoundedVector2 =             N(2)*data.m_ext(2);
const double crRightHandSideBoundedVector3 =             N(3)*data.m_ext(3);
const double crRightHandSideBoundedVector4 =             crRightHandSideBoundedVector0 + crRightHandSideBoundedVector1 + crRightHandSideBoundedVector2 + crRightHandSideBoundedVector3;
const double crRightHandSideBoundedVector5 =             N(0)*crRightHandSideBoundedVector4;
const double crRightHandSideBoundedVector6 =             DN_DX(0,1)*data.U(0,2);
const double crRightHandSideBoundedVector7 =             DN_DX(1,1)*data.U(1,2);
const double crRightHandSideBoundedVector8 =             DN_DX(2,1)*data.U(2,2);
const double crRightHandSideBoundedVector9 =             DN_DX(3,1)*data.U(3,2);
const double crRightHandSideBoundedVector10 =             crRightHandSideBoundedVector6 + crRightHandSideBoundedVector7 + crRightHandSideBoundedVector8 + crRightHandSideBoundedVector9;
const double crRightHandSideBoundedVector11 =             DN_DX(0,0)*data.U(0,1);
const double crRightHandSideBoundedVector12 =             DN_DX(1,0)*data.U(1,1);
const double crRightHandSideBoundedVector13 =             DN_DX(2,0)*data.U(2,1);
const double crRightHandSideBoundedVector14 =             DN_DX(3,0)*data.U(3,1);
const double crRightHandSideBoundedVector15 =             crRightHandSideBoundedVector11 + crRightHandSideBoundedVector12 + crRightHandSideBoundedVector13 + crRightHandSideBoundedVector14;
const double crRightHandSideBoundedVector16 =             crRightHandSideBoundedVector10 + crRightHandSideBoundedVector15;
const double crRightHandSideBoundedVector17 =             DN_DX(0,0)*data.U(0,0) + DN_DX(1,0)*data.U(1,0) + DN_DX(2,0)*data.U(2,0) + DN_DX(3,0)*data.U(3,0);
const double crRightHandSideBoundedVector18 =             N(0)*data.alpha_sc_nodes(0) + N(1)*data.alpha_sc_nodes(1) + N(2)*data.alpha_sc_nodes(2) + N(3)*data.alpha_sc_nodes(3);
const double crRightHandSideBoundedVector19 =             crRightHandSideBoundedVector17*crRightHandSideBoundedVector18;
const double crRightHandSideBoundedVector20 =             DN_DX(0,1)*data.U(0,0) + DN_DX(1,1)*data.U(1,0) + DN_DX(2,1)*data.U(2,0) + DN_DX(3,1)*data.U(3,0);
const double crRightHandSideBoundedVector21 =             crRightHandSideBoundedVector18*crRightHandSideBoundedVector20;
const double crRightHandSideBoundedVector22 =             N(0)*data.U(0,0) + N(1)*data.U(1,0) + N(2)*data.U(2,0) + N(3)*data.U(3,0);
const double crRightHandSideBoundedVector23 =             1.0/crRightHandSideBoundedVector22;
const double crRightHandSideBoundedVector24 =             fabs(crRightHandSideBoundedVector22);
const double crRightHandSideBoundedVector25 =             N(0)*data.U(0,1) + N(1)*data.U(1,1) + N(2)*data.U(2,1) + N(3)*data.U(3,1);
const double crRightHandSideBoundedVector26 =             pow(crRightHandSideBoundedVector25, 2);
const double crRightHandSideBoundedVector27 =             N(0)*data.U(0,2) + N(1)*data.U(1,2) + N(2)*data.U(2,2) + N(3)*data.U(3,2);
const double crRightHandSideBoundedVector28 =             pow(crRightHandSideBoundedVector27, 2);
const double crRightHandSideBoundedVector29 =             crRightHandSideBoundedVector26 + crRightHandSideBoundedVector28;
const double crRightHandSideBoundedVector30 =             0.5*crRightHandSideBoundedVector26;
const double crRightHandSideBoundedVector31 =             0.5*crRightHandSideBoundedVector28;
const double crRightHandSideBoundedVector32 =             N(0)*data.U(0,3);
const double crRightHandSideBoundedVector33 =             N(1)*data.U(1,3);
const double crRightHandSideBoundedVector34 =             N(2)*data.U(2,3);
const double crRightHandSideBoundedVector35 =             N(3)*data.U(3,3);
const double crRightHandSideBoundedVector36 =             -crRightHandSideBoundedVector32 - crRightHandSideBoundedVector33 - crRightHandSideBoundedVector34 - crRightHandSideBoundedVector35;
const double crRightHandSideBoundedVector37 =             crRightHandSideBoundedVector23*crRightHandSideBoundedVector30 + crRightHandSideBoundedVector23*crRightHandSideBoundedVector31 + crRightHandSideBoundedVector36;
const double crRightHandSideBoundedVector38 =             data.gamma - 1;
const double crRightHandSideBoundedVector39 =             crRightHandSideBoundedVector23*crRightHandSideBoundedVector38;
const double crRightHandSideBoundedVector40 =             crRightHandSideBoundedVector37*crRightHandSideBoundedVector39;
const double crRightHandSideBoundedVector41 =             1.0/data.gamma;
const double crRightHandSideBoundedVector42 =             N(0)*data.r_ext(0) + N(1)*data.r_ext(1) + N(2)*data.r_ext(2) + N(3)*data.r_ext(3);
const double crRightHandSideBoundedVector43 =             pow(crRightHandSideBoundedVector42, 2);
const double crRightHandSideBoundedVector44 =             N(0)*data.f_ext(0,0) + N(1)*data.f_ext(1,0) + N(2)*data.f_ext(2,0) + N(3)*data.f_ext(3,0);
const double crRightHandSideBoundedVector45 =             N(0)*data.f_ext(0,1) + N(1)*data.f_ext(1,1) + N(2)*data.f_ext(2,1) + N(3)*data.f_ext(3,1);
const double crRightHandSideBoundedVector46 =             crRightHandSideBoundedVector40*data.gamma*(pow(crRightHandSideBoundedVector44, 2) + pow(crRightHandSideBoundedVector45, 2));
const double crRightHandSideBoundedVector47 =             0.70710678118654757*crRightHandSideBoundedVector24*crRightHandSideBoundedVector41*stab_c3*sqrt((crRightHandSideBoundedVector43 - 2.0*crRightHandSideBoundedVector46 + 2.0*sqrt(0.25*crRightHandSideBoundedVector43 - crRightHandSideBoundedVector46)*fabs(crRightHandSideBoundedVector42))/(pow(crRightHandSideBoundedVector37, 2)*pow(crRightHandSideBoundedVector38, 2))) + stab_c2*(sqrt(data.gamma)*sqrt(-crRightHandSideBoundedVector40) + sqrt(crRightHandSideBoundedVector29)/crRightHandSideBoundedVector24)/data.h;
const double crRightHandSideBoundedVector48 =             1.0*(N(0)*data.ResProj(0,0) + N(0)*data.dUdt(0,0) + N(1)*data.ResProj(1,0) + N(1)*data.dUdt(1,0) + N(2)*data.ResProj(2,0) + N(2)*data.dUdt(2,0) + N(3)*data.ResProj(3,0) + N(3)*data.dUdt(3,0) - crRightHandSideBoundedVector0 - crRightHandSideBoundedVector1 + crRightHandSideBoundedVector16 - crRightHandSideBoundedVector2 - crRightHandSideBoundedVector3)/crRightHandSideBoundedVector47;
const double crRightHandSideBoundedVector49 =             crRightHandSideBoundedVector23*crRightHandSideBoundedVector48;
const double crRightHandSideBoundedVector50 =             crRightHandSideBoundedVector22*crRightHandSideBoundedVector44;
const double crRightHandSideBoundedVector51 =             crRightHandSideBoundedVector15*crRightHandSideBoundedVector25;
const double crRightHandSideBoundedVector52 =             1.0*data.gamma;
const double crRightHandSideBoundedVector53 =             crRightHandSideBoundedVector23*(3.0 - crRightHandSideBoundedVector52);
const double crRightHandSideBoundedVector54 =             DN_DX(0,0)*data.U(0,2);
const double crRightHandSideBoundedVector55 =             DN_DX(1,0)*data.U(1,2);
const double crRightHandSideBoundedVector56 =             DN_DX(2,0)*data.U(2,2);
const double crRightHandSideBoundedVector57 =             DN_DX(3,0)*data.U(3,2);
const double crRightHandSideBoundedVector58 =             crRightHandSideBoundedVector54 + crRightHandSideBoundedVector55 + crRightHandSideBoundedVector56 + crRightHandSideBoundedVector57;
const double crRightHandSideBoundedVector59 =             crRightHandSideBoundedVector27*crRightHandSideBoundedVector58;
const double crRightHandSideBoundedVector60 =             1.0*crRightHandSideBoundedVector39;
const double crRightHandSideBoundedVector61 =             -crRightHandSideBoundedVector59*crRightHandSideBoundedVector60;
const double crRightHandSideBoundedVector62 =             DN_DX(0,0)*data.U(0,3) + DN_DX(1,0)*data.U(1,3) + DN_DX(2,0)*data.U(2,3) + DN_DX(3,0)*data.U(3,3);
const double crRightHandSideBoundedVector63 =             DN_DX(0,1)*data.U(0,1);
const double crRightHandSideBoundedVector64 =             DN_DX(1,1)*data.U(1,1);
const double crRightHandSideBoundedVector65 =             DN_DX(2,1)*data.U(2,1);
const double crRightHandSideBoundedVector66 =             DN_DX(3,1)*data.U(3,1);
const double crRightHandSideBoundedVector67 =             crRightHandSideBoundedVector63 + crRightHandSideBoundedVector64 + crRightHandSideBoundedVector65 + crRightHandSideBoundedVector66;
const double crRightHandSideBoundedVector68 =             crRightHandSideBoundedVector27*crRightHandSideBoundedVector67;
const double crRightHandSideBoundedVector69 =             crRightHandSideBoundedVector10*crRightHandSideBoundedVector25;
const double crRightHandSideBoundedVector70 =             pow(crRightHandSideBoundedVector22, -2);
const double crRightHandSideBoundedVector71 =             0.5*crRightHandSideBoundedVector29;
const double crRightHandSideBoundedVector72 =             crRightHandSideBoundedVector38*crRightHandSideBoundedVector71;
const double crRightHandSideBoundedVector73 =             -crRightHandSideBoundedVector26 + crRightHandSideBoundedVector72;
const double crRightHandSideBoundedVector74 =             crRightHandSideBoundedVector70*crRightHandSideBoundedVector73;
const double crRightHandSideBoundedVector75 =             crRightHandSideBoundedVector20*crRightHandSideBoundedVector27;
const double crRightHandSideBoundedVector76 =             crRightHandSideBoundedVector70*crRightHandSideBoundedVector75;
const double crRightHandSideBoundedVector77 =             crRightHandSideBoundedVector17*crRightHandSideBoundedVector74 + crRightHandSideBoundedVector23*crRightHandSideBoundedVector68 + crRightHandSideBoundedVector23*crRightHandSideBoundedVector69 - crRightHandSideBoundedVector25*crRightHandSideBoundedVector76 + crRightHandSideBoundedVector38*crRightHandSideBoundedVector62 + crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector78 =             N(0)*data.mu_sc_nodes(0);
const double crRightHandSideBoundedVector79 =             N(1)*data.mu_sc_nodes(1);
const double crRightHandSideBoundedVector80 =             N(2)*data.mu_sc_nodes(2);
const double crRightHandSideBoundedVector81 =             N(3)*data.mu_sc_nodes(3);
const double crRightHandSideBoundedVector82 =             crRightHandSideBoundedVector78 + crRightHandSideBoundedVector79 + crRightHandSideBoundedVector80 + crRightHandSideBoundedVector81 + data.mu;
const double crRightHandSideBoundedVector83 =             crRightHandSideBoundedVector23*stab_c1/pow(data.h, 2);
const double crRightHandSideBoundedVector84 =             1.0/(crRightHandSideBoundedVector47 + 1.3333333333333333*crRightHandSideBoundedVector82*crRightHandSideBoundedVector83);
const double crRightHandSideBoundedVector85 =             crRightHandSideBoundedVector84*(N(0)*data.ResProj(0,1) + N(0)*data.dUdt(0,1) + N(1)*data.ResProj(1,1) + N(1)*data.dUdt(1,1) + N(2)*data.ResProj(2,1) + N(2)*data.dUdt(2,1) + N(3)*data.ResProj(3,1) + N(3)*data.dUdt(3,1) - crRightHandSideBoundedVector50 + crRightHandSideBoundedVector51*crRightHandSideBoundedVector53 + crRightHandSideBoundedVector77);
const double crRightHandSideBoundedVector86 =             crRightHandSideBoundedVector22*crRightHandSideBoundedVector45;
const double crRightHandSideBoundedVector87 =             crRightHandSideBoundedVector10*crRightHandSideBoundedVector27;
const double crRightHandSideBoundedVector88 =             crRightHandSideBoundedVector25*crRightHandSideBoundedVector67;
const double crRightHandSideBoundedVector89 =             -crRightHandSideBoundedVector60*crRightHandSideBoundedVector88;
const double crRightHandSideBoundedVector90 =             DN_DX(0,1)*data.U(0,3) + DN_DX(1,1)*data.U(1,3) + DN_DX(2,1)*data.U(2,3) + DN_DX(3,1)*data.U(3,3);
const double crRightHandSideBoundedVector91 =             crRightHandSideBoundedVector15*crRightHandSideBoundedVector27;
const double crRightHandSideBoundedVector92 =             crRightHandSideBoundedVector25*crRightHandSideBoundedVector58;
const double crRightHandSideBoundedVector93 =             -crRightHandSideBoundedVector28 + crRightHandSideBoundedVector72;
const double crRightHandSideBoundedVector94 =             crRightHandSideBoundedVector70*crRightHandSideBoundedVector93;
const double crRightHandSideBoundedVector95 =             crRightHandSideBoundedVector17*crRightHandSideBoundedVector25;
const double crRightHandSideBoundedVector96 =             crRightHandSideBoundedVector27*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector97 =             crRightHandSideBoundedVector20*crRightHandSideBoundedVector94 + crRightHandSideBoundedVector23*crRightHandSideBoundedVector91 + crRightHandSideBoundedVector23*crRightHandSideBoundedVector92 + crRightHandSideBoundedVector38*crRightHandSideBoundedVector90 + crRightHandSideBoundedVector89 - crRightHandSideBoundedVector95*crRightHandSideBoundedVector96;
const double crRightHandSideBoundedVector98 =             crRightHandSideBoundedVector84*(N(0)*data.ResProj(0,2) + N(0)*data.dUdt(0,2) + N(1)*data.ResProj(1,2) + N(1)*data.dUdt(1,2) + N(2)*data.ResProj(2,2) + N(2)*data.dUdt(2,2) + N(3)*data.ResProj(3,2) + N(3)*data.dUdt(3,2) + crRightHandSideBoundedVector53*crRightHandSideBoundedVector87 - crRightHandSideBoundedVector86 + crRightHandSideBoundedVector97);
const double crRightHandSideBoundedVector99 =             DN_DX(0,1)*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector100 =             crRightHandSideBoundedVector17*crRightHandSideBoundedVector27;
const double crRightHandSideBoundedVector101 =             crRightHandSideBoundedVector20*crRightHandSideBoundedVector25;
const double crRightHandSideBoundedVector102 =             crRightHandSideBoundedVector82*(-crRightHandSideBoundedVector100 - crRightHandSideBoundedVector101 + crRightHandSideBoundedVector22*crRightHandSideBoundedVector58 + crRightHandSideBoundedVector22*crRightHandSideBoundedVector67);
const double crRightHandSideBoundedVector103 =             crRightHandSideBoundedVector15*crRightHandSideBoundedVector22 - crRightHandSideBoundedVector95;
const double crRightHandSideBoundedVector104 =             2*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector105 =             crRightHandSideBoundedVector10*crRightHandSideBoundedVector22 - crRightHandSideBoundedVector75;
const double crRightHandSideBoundedVector106 =             (crRightHandSideBoundedVector103 + crRightHandSideBoundedVector105)*(N(0)*data.beta_sc_nodes(0) + N(1)*data.beta_sc_nodes(1) + N(2)*data.beta_sc_nodes(2) + N(3)*data.beta_sc_nodes(3) - 0.66666666666666663*crRightHandSideBoundedVector78 - 0.66666666666666663*crRightHandSideBoundedVector79 - 0.66666666666666663*crRightHandSideBoundedVector80 - 0.66666666666666663*crRightHandSideBoundedVector81 - 0.66666666666666663*data.mu);
const double crRightHandSideBoundedVector107 =             crRightHandSideBoundedVector103*crRightHandSideBoundedVector104 + crRightHandSideBoundedVector106;
const double crRightHandSideBoundedVector108 =             DN_DX(0,0)*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector109 =             crRightHandSideBoundedVector52 - 3.0;
const double crRightHandSideBoundedVector110 =             crRightHandSideBoundedVector109*crRightHandSideBoundedVector51;
const double crRightHandSideBoundedVector111 =             -crRightHandSideBoundedVector110*crRightHandSideBoundedVector23 + crRightHandSideBoundedVector77;
const double crRightHandSideBoundedVector112 =             N(0)*crRightHandSideBoundedVector44;
const double crRightHandSideBoundedVector113 =             DN_DX(0,1)*crRightHandSideBoundedVector25;
const double crRightHandSideBoundedVector114 =             crRightHandSideBoundedVector113*crRightHandSideBoundedVector96;
const double crRightHandSideBoundedVector115 =             crRightHandSideBoundedVector23*crRightHandSideBoundedVector75;
const double crRightHandSideBoundedVector116 =             2*crRightHandSideBoundedVector115;
const double crRightHandSideBoundedVector117 =             -crRightHandSideBoundedVector116*crRightHandSideBoundedVector25 + crRightHandSideBoundedVector68 + crRightHandSideBoundedVector69;
const double crRightHandSideBoundedVector118 =             N(0)*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector119 =             1.0*crRightHandSideBoundedVector38;
const double crRightHandSideBoundedVector120 =             crRightHandSideBoundedVector17*crRightHandSideBoundedVector23;
const double crRightHandSideBoundedVector121 =             crRightHandSideBoundedVector110 + crRightHandSideBoundedVector119*crRightHandSideBoundedVector59 - 2*crRightHandSideBoundedVector120*crRightHandSideBoundedVector73;
const double crRightHandSideBoundedVector122 =             DN_DX(0,0)*crRightHandSideBoundedVector27;
const double crRightHandSideBoundedVector123 =             crRightHandSideBoundedVector101*crRightHandSideBoundedVector23 - crRightHandSideBoundedVector63 - crRightHandSideBoundedVector64 - crRightHandSideBoundedVector65 - crRightHandSideBoundedVector66;
const double crRightHandSideBoundedVector124 =             N(0)*crRightHandSideBoundedVector123;
const double crRightHandSideBoundedVector125 =             crRightHandSideBoundedVector52 - 1.0;
const double crRightHandSideBoundedVector126 =             crRightHandSideBoundedVector100*crRightHandSideBoundedVector23 - crRightHandSideBoundedVector54 - crRightHandSideBoundedVector55 - crRightHandSideBoundedVector56 - crRightHandSideBoundedVector57;
const double crRightHandSideBoundedVector127 =             N(0)*crRightHandSideBoundedVector126;
const double crRightHandSideBoundedVector128 =             crRightHandSideBoundedVector23*crRightHandSideBoundedVector98;
const double crRightHandSideBoundedVector129 =             crRightHandSideBoundedVector23*crRightHandSideBoundedVector95;
const double crRightHandSideBoundedVector130 =             -crRightHandSideBoundedVector11 - crRightHandSideBoundedVector12 + crRightHandSideBoundedVector129 - crRightHandSideBoundedVector13 - crRightHandSideBoundedVector14;
const double crRightHandSideBoundedVector131 =             N(0)*crRightHandSideBoundedVector130;
const double crRightHandSideBoundedVector132 =             DN_DX(0,0)*crRightHandSideBoundedVector25;
const double crRightHandSideBoundedVector133 =             DN_DX(0,1)*crRightHandSideBoundedVector27;
const double crRightHandSideBoundedVector134 =             crRightHandSideBoundedVector115 - crRightHandSideBoundedVector6 - crRightHandSideBoundedVector7 - crRightHandSideBoundedVector8 - crRightHandSideBoundedVector9;
const double crRightHandSideBoundedVector135 =             N(0)*crRightHandSideBoundedVector134;
const double crRightHandSideBoundedVector136 =             crRightHandSideBoundedVector133 - crRightHandSideBoundedVector135;
const double crRightHandSideBoundedVector137 =             crRightHandSideBoundedVector23*crRightHandSideBoundedVector85;
const double crRightHandSideBoundedVector138 =             (N(0)*data.lamb_sc_nodes(0) + N(1)*data.lamb_sc_nodes(1) + N(2)*data.lamb_sc_nodes(2) + N(3)*data.lamb_sc_nodes(3) + data.lambda)/data.c_v;
const double crRightHandSideBoundedVector139 =             crRightHandSideBoundedVector22*crRightHandSideBoundedVector42;
const double crRightHandSideBoundedVector140 =             crRightHandSideBoundedVector25*crRightHandSideBoundedVector44;
const double crRightHandSideBoundedVector141 =             crRightHandSideBoundedVector27*crRightHandSideBoundedVector45;
const double crRightHandSideBoundedVector142 =             crRightHandSideBoundedVector62*data.gamma;
const double crRightHandSideBoundedVector143 =             crRightHandSideBoundedVector23*crRightHandSideBoundedVector25;
const double crRightHandSideBoundedVector144 =             crRightHandSideBoundedVector90*data.gamma;
const double crRightHandSideBoundedVector145 =             crRightHandSideBoundedVector23*crRightHandSideBoundedVector27;
const double crRightHandSideBoundedVector146 =             1.0*crRightHandSideBoundedVector92;
const double crRightHandSideBoundedVector147 =             1.0*crRightHandSideBoundedVector68;
const double crRightHandSideBoundedVector148 =             crRightHandSideBoundedVector147*crRightHandSideBoundedVector25;
const double crRightHandSideBoundedVector149 =             crRightHandSideBoundedVector26*crRightHandSideBoundedVector60;
const double crRightHandSideBoundedVector150 =             crRightHandSideBoundedVector32 + crRightHandSideBoundedVector33 + crRightHandSideBoundedVector34 + crRightHandSideBoundedVector35;
const double crRightHandSideBoundedVector151 =             crRightHandSideBoundedVector38*(crRightHandSideBoundedVector150 - crRightHandSideBoundedVector23*(crRightHandSideBoundedVector30 + crRightHandSideBoundedVector31));
const double crRightHandSideBoundedVector152 =             crRightHandSideBoundedVector150 + crRightHandSideBoundedVector151;
const double crRightHandSideBoundedVector153 =             crRightHandSideBoundedVector28*crRightHandSideBoundedVector60;
const double crRightHandSideBoundedVector154 =             -crRightHandSideBoundedVector151 + crRightHandSideBoundedVector36;
const double crRightHandSideBoundedVector155 =             crRightHandSideBoundedVector154 + crRightHandSideBoundedVector39*crRightHandSideBoundedVector71;
const double crRightHandSideBoundedVector156 =             crRightHandSideBoundedVector155*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector157 =             (N(0)*data.ResProj(0,3) + N(0)*data.dUdt(0,3) + N(1)*data.ResProj(1,3) + N(1)*data.dUdt(1,3) + N(2)*data.ResProj(2,3) + N(2)*data.dUdt(2,3) + N(3)*data.ResProj(3,3) + N(3)*data.dUdt(3,3) + crRightHandSideBoundedVector10*crRightHandSideBoundedVector23*(crRightHandSideBoundedVector152 - crRightHandSideBoundedVector153) - crRightHandSideBoundedVector139 - crRightHandSideBoundedVector140 - crRightHandSideBoundedVector141 + crRightHandSideBoundedVector142*crRightHandSideBoundedVector143 + crRightHandSideBoundedVector144*crRightHandSideBoundedVector145 - crRightHandSideBoundedVector146*crRightHandSideBoundedVector38*crRightHandSideBoundedVector96 - crRightHandSideBoundedVector148*crRightHandSideBoundedVector38*crRightHandSideBoundedVector70 + crRightHandSideBoundedVector15*crRightHandSideBoundedVector23*(-crRightHandSideBoundedVector149 + crRightHandSideBoundedVector152) + crRightHandSideBoundedVector155*crRightHandSideBoundedVector76 + crRightHandSideBoundedVector156*crRightHandSideBoundedVector95)/(crRightHandSideBoundedVector138*crRightHandSideBoundedVector41*crRightHandSideBoundedVector83 + crRightHandSideBoundedVector47);
const double crRightHandSideBoundedVector158 =             crRightHandSideBoundedVector119*crRightHandSideBoundedVector157;
const double crRightHandSideBoundedVector159 =             crRightHandSideBoundedVector104*crRightHandSideBoundedVector105 + crRightHandSideBoundedVector106;
const double crRightHandSideBoundedVector160 =             crRightHandSideBoundedVector109*crRightHandSideBoundedVector87;
const double crRightHandSideBoundedVector161 =             -crRightHandSideBoundedVector160*crRightHandSideBoundedVector23 + crRightHandSideBoundedVector97;
const double crRightHandSideBoundedVector162 =             N(0)*crRightHandSideBoundedVector45;
const double crRightHandSideBoundedVector163 =             crRightHandSideBoundedVector132*crRightHandSideBoundedVector96;
const double crRightHandSideBoundedVector164 =             2*crRightHandSideBoundedVector129;
const double crRightHandSideBoundedVector165 =             -crRightHandSideBoundedVector164*crRightHandSideBoundedVector27 + crRightHandSideBoundedVector91 + crRightHandSideBoundedVector92;
const double crRightHandSideBoundedVector166 =             crRightHandSideBoundedVector20*crRightHandSideBoundedVector23;
const double crRightHandSideBoundedVector167 =             crRightHandSideBoundedVector119*crRightHandSideBoundedVector88 + crRightHandSideBoundedVector160 - 2*crRightHandSideBoundedVector166*crRightHandSideBoundedVector93;
const double crRightHandSideBoundedVector168 =             -crRightHandSideBoundedVector131 + crRightHandSideBoundedVector132;
const double crRightHandSideBoundedVector169 =             crRightHandSideBoundedVector139 + crRightHandSideBoundedVector140 + crRightHandSideBoundedVector141;
const double crRightHandSideBoundedVector170 =             crRightHandSideBoundedVector102*crRightHandSideBoundedVector145 + crRightHandSideBoundedVector107*crRightHandSideBoundedVector143 + crRightHandSideBoundedVector138*(crRightHandSideBoundedVector120*crRightHandSideBoundedVector26 + crRightHandSideBoundedVector120*crRightHandSideBoundedVector28 - crRightHandSideBoundedVector150*crRightHandSideBoundedVector17 + crRightHandSideBoundedVector22*crRightHandSideBoundedVector62 - crRightHandSideBoundedVector51 - crRightHandSideBoundedVector59);
const double crRightHandSideBoundedVector171 =             crRightHandSideBoundedVector102*crRightHandSideBoundedVector143 + crRightHandSideBoundedVector138*(-crRightHandSideBoundedVector150*crRightHandSideBoundedVector20 + crRightHandSideBoundedVector166*crRightHandSideBoundedVector26 + crRightHandSideBoundedVector166*crRightHandSideBoundedVector28 + crRightHandSideBoundedVector22*crRightHandSideBoundedVector90 - crRightHandSideBoundedVector87 - crRightHandSideBoundedVector88) + crRightHandSideBoundedVector145*crRightHandSideBoundedVector159;
const double crRightHandSideBoundedVector172 =             crRightHandSideBoundedVector142*crRightHandSideBoundedVector25;
const double crRightHandSideBoundedVector173 =             crRightHandSideBoundedVector144*crRightHandSideBoundedVector27;
const double crRightHandSideBoundedVector174 =             crRightHandSideBoundedVector149 + crRightHandSideBoundedVector154;
const double crRightHandSideBoundedVector175 =             crRightHandSideBoundedVector153 + crRightHandSideBoundedVector154;
const double crRightHandSideBoundedVector176 =             -crRightHandSideBoundedVector10*crRightHandSideBoundedVector175 + crRightHandSideBoundedVector115*crRightHandSideBoundedVector155 + crRightHandSideBoundedVector129*crRightHandSideBoundedVector155 - crRightHandSideBoundedVector146*crRightHandSideBoundedVector27*crRightHandSideBoundedVector39 - crRightHandSideBoundedVector148*crRightHandSideBoundedVector39 - crRightHandSideBoundedVector15*crRightHandSideBoundedVector174 + crRightHandSideBoundedVector172 + crRightHandSideBoundedVector173;
const double crRightHandSideBoundedVector177 =             N(0)*crRightHandSideBoundedVector23;
const double crRightHandSideBoundedVector178 =             -2.0*crRightHandSideBoundedVector115*crRightHandSideBoundedVector25 + crRightHandSideBoundedVector147 + 1.0*crRightHandSideBoundedVector69;
const double crRightHandSideBoundedVector179 =             crRightHandSideBoundedVector118*crRightHandSideBoundedVector38;
const double crRightHandSideBoundedVector180 =             crRightHandSideBoundedVector174*crRightHandSideBoundedVector23;
const double crRightHandSideBoundedVector181 =             3.0*crRightHandSideBoundedVector39;
const double crRightHandSideBoundedVector182 =             2.0*crRightHandSideBoundedVector39;
const double crRightHandSideBoundedVector183 =             crRightHandSideBoundedVector155 + crRightHandSideBoundedVector182*crRightHandSideBoundedVector26;
const double crRightHandSideBoundedVector184 =             crRightHandSideBoundedVector120*crRightHandSideBoundedVector183 + crRightHandSideBoundedVector142 - crRightHandSideBoundedVector181*crRightHandSideBoundedVector51 + crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector185 =             -2.0*crRightHandSideBoundedVector129*crRightHandSideBoundedVector27 + crRightHandSideBoundedVector146 + 1.0*crRightHandSideBoundedVector91;
const double crRightHandSideBoundedVector186 =             crRightHandSideBoundedVector175*crRightHandSideBoundedVector23;
const double crRightHandSideBoundedVector187 =             crRightHandSideBoundedVector155 + crRightHandSideBoundedVector182*crRightHandSideBoundedVector28;
const double crRightHandSideBoundedVector188 =             crRightHandSideBoundedVector144 + crRightHandSideBoundedVector166*crRightHandSideBoundedVector187 - crRightHandSideBoundedVector181*crRightHandSideBoundedVector87 + crRightHandSideBoundedVector89;
const double crRightHandSideBoundedVector189 =             crRightHandSideBoundedVector157*crRightHandSideBoundedVector23*crRightHandSideBoundedVector52;
const double crRightHandSideBoundedVector190 =             crRightHandSideBoundedVector154 + crRightHandSideBoundedVector29*crRightHandSideBoundedVector39;
const double crRightHandSideBoundedVector191 =             crRightHandSideBoundedVector15*crRightHandSideBoundedVector183 - crRightHandSideBoundedVector164*crRightHandSideBoundedVector190 - crRightHandSideBoundedVector172 + crRightHandSideBoundedVector182*crRightHandSideBoundedVector27*crRightHandSideBoundedVector92;
const double crRightHandSideBoundedVector192 =             crRightHandSideBoundedVector10*crRightHandSideBoundedVector187 - crRightHandSideBoundedVector116*crRightHandSideBoundedVector190 - crRightHandSideBoundedVector173 + crRightHandSideBoundedVector182*crRightHandSideBoundedVector25*crRightHandSideBoundedVector68;
const double crRightHandSideBoundedVector193 =             N(1)*crRightHandSideBoundedVector4;
const double crRightHandSideBoundedVector194 =             DN_DX(1,1)*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector195 =             DN_DX(1,0)*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector196 =             N(1)*crRightHandSideBoundedVector44;
const double crRightHandSideBoundedVector197 =             DN_DX(1,1)*crRightHandSideBoundedVector25;
const double crRightHandSideBoundedVector198 =             crRightHandSideBoundedVector197*crRightHandSideBoundedVector96;
const double crRightHandSideBoundedVector199 =             N(1)*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector200 =             DN_DX(1,0)*crRightHandSideBoundedVector27;
const double crRightHandSideBoundedVector201 =             N(1)*crRightHandSideBoundedVector123;
const double crRightHandSideBoundedVector202 =             N(1)*crRightHandSideBoundedVector126;
const double crRightHandSideBoundedVector203 =             N(1)*crRightHandSideBoundedVector130;
const double crRightHandSideBoundedVector204 =             DN_DX(1,0)*crRightHandSideBoundedVector25;
const double crRightHandSideBoundedVector205 =             DN_DX(1,1)*crRightHandSideBoundedVector27;
const double crRightHandSideBoundedVector206 =             N(1)*crRightHandSideBoundedVector134;
const double crRightHandSideBoundedVector207 =             crRightHandSideBoundedVector205 - crRightHandSideBoundedVector206;
const double crRightHandSideBoundedVector208 =             N(1)*crRightHandSideBoundedVector45;
const double crRightHandSideBoundedVector209 =             crRightHandSideBoundedVector204*crRightHandSideBoundedVector96;
const double crRightHandSideBoundedVector210 =             -crRightHandSideBoundedVector203 + crRightHandSideBoundedVector204;
const double crRightHandSideBoundedVector211 =             N(1)*crRightHandSideBoundedVector23;
const double crRightHandSideBoundedVector212 =             crRightHandSideBoundedVector199*crRightHandSideBoundedVector38;
const double crRightHandSideBoundedVector213 =             N(2)*crRightHandSideBoundedVector4;
const double crRightHandSideBoundedVector214 =             DN_DX(2,1)*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector215 =             DN_DX(2,0)*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector216 =             N(2)*crRightHandSideBoundedVector44;
const double crRightHandSideBoundedVector217 =             DN_DX(2,1)*crRightHandSideBoundedVector25;
const double crRightHandSideBoundedVector218 =             crRightHandSideBoundedVector217*crRightHandSideBoundedVector96;
const double crRightHandSideBoundedVector219 =             N(2)*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector220 =             DN_DX(2,0)*crRightHandSideBoundedVector27;
const double crRightHandSideBoundedVector221 =             N(2)*crRightHandSideBoundedVector123;
const double crRightHandSideBoundedVector222 =             N(2)*crRightHandSideBoundedVector126;
const double crRightHandSideBoundedVector223 =             N(2)*crRightHandSideBoundedVector130;
const double crRightHandSideBoundedVector224 =             DN_DX(2,0)*crRightHandSideBoundedVector25;
const double crRightHandSideBoundedVector225 =             DN_DX(2,1)*crRightHandSideBoundedVector27;
const double crRightHandSideBoundedVector226 =             N(2)*crRightHandSideBoundedVector134;
const double crRightHandSideBoundedVector227 =             crRightHandSideBoundedVector225 - crRightHandSideBoundedVector226;
const double crRightHandSideBoundedVector228 =             N(2)*crRightHandSideBoundedVector45;
const double crRightHandSideBoundedVector229 =             crRightHandSideBoundedVector224*crRightHandSideBoundedVector96;
const double crRightHandSideBoundedVector230 =             -crRightHandSideBoundedVector223 + crRightHandSideBoundedVector224;
const double crRightHandSideBoundedVector231 =             N(2)*crRightHandSideBoundedVector23;
const double crRightHandSideBoundedVector232 =             crRightHandSideBoundedVector219*crRightHandSideBoundedVector38;
const double crRightHandSideBoundedVector233 =             N(3)*crRightHandSideBoundedVector4;
const double crRightHandSideBoundedVector234 =             DN_DX(3,1)*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector235 =             DN_DX(3,0)*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector236 =             N(3)*crRightHandSideBoundedVector44;
const double crRightHandSideBoundedVector237 =             DN_DX(3,1)*crRightHandSideBoundedVector25;
const double crRightHandSideBoundedVector238 =             crRightHandSideBoundedVector237*crRightHandSideBoundedVector96;
const double crRightHandSideBoundedVector239 =             N(3)*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector240 =             DN_DX(3,0)*crRightHandSideBoundedVector27;
const double crRightHandSideBoundedVector241 =             N(3)*crRightHandSideBoundedVector123;
const double crRightHandSideBoundedVector242 =             N(3)*crRightHandSideBoundedVector126;
const double crRightHandSideBoundedVector243 =             N(3)*crRightHandSideBoundedVector130;
const double crRightHandSideBoundedVector244 =             DN_DX(3,0)*crRightHandSideBoundedVector25;
const double crRightHandSideBoundedVector245 =             DN_DX(3,1)*crRightHandSideBoundedVector27;
const double crRightHandSideBoundedVector246 =             N(3)*crRightHandSideBoundedVector134;
const double crRightHandSideBoundedVector247 =             crRightHandSideBoundedVector245 - crRightHandSideBoundedVector246;
const double crRightHandSideBoundedVector248 =             N(3)*crRightHandSideBoundedVector45;
const double crRightHandSideBoundedVector249 =             crRightHandSideBoundedVector244*crRightHandSideBoundedVector96;
const double crRightHandSideBoundedVector250 =             -crRightHandSideBoundedVector243 + crRightHandSideBoundedVector244;
const double crRightHandSideBoundedVector251 =             N(3)*crRightHandSideBoundedVector23;
const double crRightHandSideBoundedVector252 =             crRightHandSideBoundedVector239*crRightHandSideBoundedVector38;
            rRightHandSideBoundedVector[0] += DN_DX(0,0)*crRightHandSideBoundedVector19 - DN_DX(0,0)*crRightHandSideBoundedVector85 + DN_DX(0,1)*crRightHandSideBoundedVector21 - DN_DX(0,1)*crRightHandSideBoundedVector98 - N(0)*crRightHandSideBoundedVector16 - crRightHandSideBoundedVector49*crRightHandSideBoundedVector5 + crRightHandSideBoundedVector5;
            rRightHandSideBoundedVector[1] += -DN_DX(0,0)*crRightHandSideBoundedVector158 - N(0)*crRightHandSideBoundedVector111 + N(0)*crRightHandSideBoundedVector50 - crRightHandSideBoundedVector102*crRightHandSideBoundedVector99 - crRightHandSideBoundedVector107*crRightHandSideBoundedVector108 - crRightHandSideBoundedVector128*(crRightHandSideBoundedVector113 - crRightHandSideBoundedVector119*crRightHandSideBoundedVector122 - crRightHandSideBoundedVector124 + crRightHandSideBoundedVector125*crRightHandSideBoundedVector127) - crRightHandSideBoundedVector137*(crRightHandSideBoundedVector109*crRightHandSideBoundedVector131 - crRightHandSideBoundedVector109*crRightHandSideBoundedVector132 + crRightHandSideBoundedVector136) - crRightHandSideBoundedVector48*(DN_DX(0,0)*crRightHandSideBoundedVector74 + crRightHandSideBoundedVector112 - crRightHandSideBoundedVector114 - crRightHandSideBoundedVector117*crRightHandSideBoundedVector118 + crRightHandSideBoundedVector118*crRightHandSideBoundedVector121);
            rRightHandSideBoundedVector[2] += -DN_DX(0,1)*crRightHandSideBoundedVector158 - N(0)*crRightHandSideBoundedVector161 + N(0)*crRightHandSideBoundedVector86 - crRightHandSideBoundedVector102*crRightHandSideBoundedVector108 - crRightHandSideBoundedVector128*(-crRightHandSideBoundedVector109*crRightHandSideBoundedVector133 + crRightHandSideBoundedVector109*crRightHandSideBoundedVector135 + crRightHandSideBoundedVector168) - crRightHandSideBoundedVector137*(-crRightHandSideBoundedVector113*crRightHandSideBoundedVector119 + crRightHandSideBoundedVector122 + crRightHandSideBoundedVector124*crRightHandSideBoundedVector125 - crRightHandSideBoundedVector127) - crRightHandSideBoundedVector159*crRightHandSideBoundedVector99 - crRightHandSideBoundedVector48*(DN_DX(0,1)*crRightHandSideBoundedVector94 - crRightHandSideBoundedVector118*crRightHandSideBoundedVector165 + crRightHandSideBoundedVector118*crRightHandSideBoundedVector167 + crRightHandSideBoundedVector162 - crRightHandSideBoundedVector163);
            rRightHandSideBoundedVector[3] += N(0)*crRightHandSideBoundedVector169 - crRightHandSideBoundedVector108*crRightHandSideBoundedVector170 - crRightHandSideBoundedVector171*crRightHandSideBoundedVector99 - crRightHandSideBoundedVector176*crRightHandSideBoundedVector177 - crRightHandSideBoundedVector189*(crRightHandSideBoundedVector136 + crRightHandSideBoundedVector168) - crRightHandSideBoundedVector48*(N(0)*crRightHandSideBoundedVector42 + crRightHandSideBoundedVector118*crRightHandSideBoundedVector191 + crRightHandSideBoundedVector118*crRightHandSideBoundedVector192 + crRightHandSideBoundedVector132*crRightHandSideBoundedVector156 + crRightHandSideBoundedVector133*crRightHandSideBoundedVector156) - crRightHandSideBoundedVector85*(-DN_DX(0,0)*crRightHandSideBoundedVector180 + crRightHandSideBoundedVector112 - crRightHandSideBoundedVector114*crRightHandSideBoundedVector119 + crRightHandSideBoundedVector177*crRightHandSideBoundedVector184 - crRightHandSideBoundedVector178*crRightHandSideBoundedVector179) - crRightHandSideBoundedVector98*(-DN_DX(0,1)*crRightHandSideBoundedVector186 - crRightHandSideBoundedVector119*crRightHandSideBoundedVector163 + crRightHandSideBoundedVector162 + crRightHandSideBoundedVector177*crRightHandSideBoundedVector188 - crRightHandSideBoundedVector179*crRightHandSideBoundedVector185);
            rRightHandSideBoundedVector[4] += DN_DX(1,0)*crRightHandSideBoundedVector19 - DN_DX(1,0)*crRightHandSideBoundedVector85 + DN_DX(1,1)*crRightHandSideBoundedVector21 - DN_DX(1,1)*crRightHandSideBoundedVector98 - N(1)*crRightHandSideBoundedVector16 - crRightHandSideBoundedVector193*crRightHandSideBoundedVector49 + crRightHandSideBoundedVector193;
            rRightHandSideBoundedVector[5] += -DN_DX(1,0)*crRightHandSideBoundedVector158 - N(1)*crRightHandSideBoundedVector111 + N(1)*crRightHandSideBoundedVector50 - crRightHandSideBoundedVector102*crRightHandSideBoundedVector194 - crRightHandSideBoundedVector107*crRightHandSideBoundedVector195 - crRightHandSideBoundedVector128*(-crRightHandSideBoundedVector119*crRightHandSideBoundedVector200 + crRightHandSideBoundedVector125*crRightHandSideBoundedVector202 + crRightHandSideBoundedVector197 - crRightHandSideBoundedVector201) - crRightHandSideBoundedVector137*(crRightHandSideBoundedVector109*crRightHandSideBoundedVector203 - crRightHandSideBoundedVector109*crRightHandSideBoundedVector204 + crRightHandSideBoundedVector207) - crRightHandSideBoundedVector48*(DN_DX(1,0)*crRightHandSideBoundedVector74 - crRightHandSideBoundedVector117*crRightHandSideBoundedVector199 + crRightHandSideBoundedVector121*crRightHandSideBoundedVector199 + crRightHandSideBoundedVector196 - crRightHandSideBoundedVector198);
            rRightHandSideBoundedVector[6] += -DN_DX(1,1)*crRightHandSideBoundedVector158 - N(1)*crRightHandSideBoundedVector161 + N(1)*crRightHandSideBoundedVector86 - crRightHandSideBoundedVector102*crRightHandSideBoundedVector195 - crRightHandSideBoundedVector128*(-crRightHandSideBoundedVector109*crRightHandSideBoundedVector205 + crRightHandSideBoundedVector109*crRightHandSideBoundedVector206 + crRightHandSideBoundedVector210) - crRightHandSideBoundedVector137*(-crRightHandSideBoundedVector119*crRightHandSideBoundedVector197 + crRightHandSideBoundedVector125*crRightHandSideBoundedVector201 + crRightHandSideBoundedVector200 - crRightHandSideBoundedVector202) - crRightHandSideBoundedVector159*crRightHandSideBoundedVector194 - crRightHandSideBoundedVector48*(DN_DX(1,1)*crRightHandSideBoundedVector94 - crRightHandSideBoundedVector165*crRightHandSideBoundedVector199 + crRightHandSideBoundedVector167*crRightHandSideBoundedVector199 + crRightHandSideBoundedVector208 - crRightHandSideBoundedVector209);
            rRightHandSideBoundedVector[7] += N(1)*crRightHandSideBoundedVector169 - crRightHandSideBoundedVector170*crRightHandSideBoundedVector195 - crRightHandSideBoundedVector171*crRightHandSideBoundedVector194 - crRightHandSideBoundedVector176*crRightHandSideBoundedVector211 - crRightHandSideBoundedVector189*(crRightHandSideBoundedVector207 + crRightHandSideBoundedVector210) - crRightHandSideBoundedVector48*(N(1)*crRightHandSideBoundedVector42 + crRightHandSideBoundedVector156*crRightHandSideBoundedVector204 + crRightHandSideBoundedVector156*crRightHandSideBoundedVector205 + crRightHandSideBoundedVector191*crRightHandSideBoundedVector199 + crRightHandSideBoundedVector192*crRightHandSideBoundedVector199) - crRightHandSideBoundedVector85*(-DN_DX(1,0)*crRightHandSideBoundedVector180 - crRightHandSideBoundedVector119*crRightHandSideBoundedVector198 - crRightHandSideBoundedVector178*crRightHandSideBoundedVector212 + crRightHandSideBoundedVector184*crRightHandSideBoundedVector211 + crRightHandSideBoundedVector196) - crRightHandSideBoundedVector98*(-DN_DX(1,1)*crRightHandSideBoundedVector186 - crRightHandSideBoundedVector119*crRightHandSideBoundedVector209 - crRightHandSideBoundedVector185*crRightHandSideBoundedVector212 + crRightHandSideBoundedVector188*crRightHandSideBoundedVector211 + crRightHandSideBoundedVector208);
            rRightHandSideBoundedVector[8] += DN_DX(2,0)*crRightHandSideBoundedVector19 - DN_DX(2,0)*crRightHandSideBoundedVector85 + DN_DX(2,1)*crRightHandSideBoundedVector21 - DN_DX(2,1)*crRightHandSideBoundedVector98 - N(2)*crRightHandSideBoundedVector16 - crRightHandSideBoundedVector213*crRightHandSideBoundedVector49 + crRightHandSideBoundedVector213;
            rRightHandSideBoundedVector[9] += -DN_DX(2,0)*crRightHandSideBoundedVector158 - N(2)*crRightHandSideBoundedVector111 + N(2)*crRightHandSideBoundedVector50 - crRightHandSideBoundedVector102*crRightHandSideBoundedVector214 - crRightHandSideBoundedVector107*crRightHandSideBoundedVector215 - crRightHandSideBoundedVector128*(-crRightHandSideBoundedVector119*crRightHandSideBoundedVector220 + crRightHandSideBoundedVector125*crRightHandSideBoundedVector222 + crRightHandSideBoundedVector217 - crRightHandSideBoundedVector221) - crRightHandSideBoundedVector137*(crRightHandSideBoundedVector109*crRightHandSideBoundedVector223 - crRightHandSideBoundedVector109*crRightHandSideBoundedVector224 + crRightHandSideBoundedVector227) - crRightHandSideBoundedVector48*(DN_DX(2,0)*crRightHandSideBoundedVector74 - crRightHandSideBoundedVector117*crRightHandSideBoundedVector219 + crRightHandSideBoundedVector121*crRightHandSideBoundedVector219 + crRightHandSideBoundedVector216 - crRightHandSideBoundedVector218);
            rRightHandSideBoundedVector[10] += -DN_DX(2,1)*crRightHandSideBoundedVector158 - N(2)*crRightHandSideBoundedVector161 + N(2)*crRightHandSideBoundedVector86 - crRightHandSideBoundedVector102*crRightHandSideBoundedVector215 - crRightHandSideBoundedVector128*(-crRightHandSideBoundedVector109*crRightHandSideBoundedVector225 + crRightHandSideBoundedVector109*crRightHandSideBoundedVector226 + crRightHandSideBoundedVector230) - crRightHandSideBoundedVector137*(-crRightHandSideBoundedVector119*crRightHandSideBoundedVector217 + crRightHandSideBoundedVector125*crRightHandSideBoundedVector221 + crRightHandSideBoundedVector220 - crRightHandSideBoundedVector222) - crRightHandSideBoundedVector159*crRightHandSideBoundedVector214 - crRightHandSideBoundedVector48*(DN_DX(2,1)*crRightHandSideBoundedVector94 - crRightHandSideBoundedVector165*crRightHandSideBoundedVector219 + crRightHandSideBoundedVector167*crRightHandSideBoundedVector219 + crRightHandSideBoundedVector228 - crRightHandSideBoundedVector229);
            rRightHandSideBoundedVector[11] += N(2)*crRightHandSideBoundedVector169 - crRightHandSideBoundedVector170*crRightHandSideBoundedVector215 - crRightHandSideBoundedVector171*crRightHandSideBoundedVector214 - crRightHandSideBoundedVector176*crRightHandSideBoundedVector231 - crRightHandSideBoundedVector189*(crRightHandSideBoundedVector227 + crRightHandSideBoundedVector230) - crRightHandSideBoundedVector48*(N(2)*crRightHandSideBoundedVector42 + crRightHandSideBoundedVector156*crRightHandSideBoundedVector224 + crRightHandSideBoundedVector156*crRightHandSideBoundedVector225 + crRightHandSideBoundedVector191*crRightHandSideBoundedVector219 + crRightHandSideBoundedVector192*crRightHandSideBoundedVector219) - crRightHandSideBoundedVector85*(-DN_DX(2,0)*crRightHandSideBoundedVector180 - crRightHandSideBoundedVector119*crRightHandSideBoundedVector218 - crRightHandSideBoundedVector178*crRightHandSideBoundedVector232 + crRightHandSideBoundedVector184*crRightHandSideBoundedVector231 + crRightHandSideBoundedVector216) - crRightHandSideBoundedVector98*(-DN_DX(2,1)*crRightHandSideBoundedVector186 - crRightHandSideBoundedVector119*crRightHandSideBoundedVector229 - crRightHandSideBoundedVector185*crRightHandSideBoundedVector232 + crRightHandSideBoundedVector188*crRightHandSideBoundedVector231 + crRightHandSideBoundedVector228);
            rRightHandSideBoundedVector[12] += DN_DX(3,0)*crRightHandSideBoundedVector19 - DN_DX(3,0)*crRightHandSideBoundedVector85 + DN_DX(3,1)*crRightHandSideBoundedVector21 - DN_DX(3,1)*crRightHandSideBoundedVector98 - N(3)*crRightHandSideBoundedVector16 - crRightHandSideBoundedVector233*crRightHandSideBoundedVector49 + crRightHandSideBoundedVector233;
            rRightHandSideBoundedVector[13] += -DN_DX(3,0)*crRightHandSideBoundedVector158 - N(3)*crRightHandSideBoundedVector111 + N(3)*crRightHandSideBoundedVector50 - crRightHandSideBoundedVector102*crRightHandSideBoundedVector234 - crRightHandSideBoundedVector107*crRightHandSideBoundedVector235 - crRightHandSideBoundedVector128*(-crRightHandSideBoundedVector119*crRightHandSideBoundedVector240 + crRightHandSideBoundedVector125*crRightHandSideBoundedVector242 + crRightHandSideBoundedVector237 - crRightHandSideBoundedVector241) - crRightHandSideBoundedVector137*(crRightHandSideBoundedVector109*crRightHandSideBoundedVector243 - crRightHandSideBoundedVector109*crRightHandSideBoundedVector244 + crRightHandSideBoundedVector247) - crRightHandSideBoundedVector48*(DN_DX(3,0)*crRightHandSideBoundedVector74 - crRightHandSideBoundedVector117*crRightHandSideBoundedVector239 + crRightHandSideBoundedVector121*crRightHandSideBoundedVector239 + crRightHandSideBoundedVector236 - crRightHandSideBoundedVector238);
            rRightHandSideBoundedVector[14] += -DN_DX(3,1)*crRightHandSideBoundedVector158 - N(3)*crRightHandSideBoundedVector161 + N(3)*crRightHandSideBoundedVector86 - crRightHandSideBoundedVector102*crRightHandSideBoundedVector235 - crRightHandSideBoundedVector128*(-crRightHandSideBoundedVector109*crRightHandSideBoundedVector245 + crRightHandSideBoundedVector109*crRightHandSideBoundedVector246 + crRightHandSideBoundedVector250) - crRightHandSideBoundedVector137*(-crRightHandSideBoundedVector119*crRightHandSideBoundedVector237 + crRightHandSideBoundedVector125*crRightHandSideBoundedVector241 + crRightHandSideBoundedVector240 - crRightHandSideBoundedVector242) - crRightHandSideBoundedVector159*crRightHandSideBoundedVector234 - crRightHandSideBoundedVector48*(DN_DX(3,1)*crRightHandSideBoundedVector94 - crRightHandSideBoundedVector165*crRightHandSideBoundedVector239 + crRightHandSideBoundedVector167*crRightHandSideBoundedVector239 + crRightHandSideBoundedVector248 - crRightHandSideBoundedVector249);
            rRightHandSideBoundedVector[15] += N(3)*crRightHandSideBoundedVector169 - crRightHandSideBoundedVector170*crRightHandSideBoundedVector235 - crRightHandSideBoundedVector171*crRightHandSideBoundedVector234 - crRightHandSideBoundedVector176*crRightHandSideBoundedVector251 - crRightHandSideBoundedVector189*(crRightHandSideBoundedVector247 + crRightHandSideBoundedVector250) - crRightHandSideBoundedVector48*(N(3)*crRightHandSideBoundedVector42 + crRightHandSideBoundedVector156*crRightHandSideBoundedVector244 + crRightHandSideBoundedVector156*crRightHandSideBoundedVector245 + crRightHandSideBoundedVector191*crRightHandSideBoundedVector239 + crRightHandSideBoundedVector192*crRightHandSideBoundedVector239) - crRightHandSideBoundedVector85*(-DN_DX(3,0)*crRightHandSideBoundedVector180 - crRightHandSideBoundedVector119*crRightHandSideBoundedVector238 - crRightHandSideBoundedVector178*crRightHandSideBoundedVector252 + crRightHandSideBoundedVector184*crRightHandSideBoundedVector251 + crRightHandSideBoundedVector236) - crRightHandSideBoundedVector98*(-DN_DX(3,1)*crRightHandSideBoundedVector186 - crRightHandSideBoundedVector119*crRightHandSideBoundedVector249 - crRightHandSideBoundedVector185*crRightHandSideBoundedVector252 + crRightHandSideBoundedVector188*crRightHandSideBoundedVector251 + crRightHandSideBoundedVector248);

        }
    }
    else
    {
        for(const auto& gauss_point: gauss_points)
        {
            r_geometry.ShapeFunctionsValues(N, gauss_point.Coordinates());
            r_geometry.InverseOfJacobian(Jinv, gauss_point.Coordinates());
            r_geometry.ShapeFunctionsLocalGradients(DN_DX_iso, gauss_point.Coordinates());
            GeometryUtils::ShapeFunctionsGradients(DN_DX_iso, Jinv, DN_DX);

const double crRightHandSideBoundedVector0 =             N(0)*data.m_ext(0);
const double crRightHandSideBoundedVector1 =             N(1)*data.m_ext(1);
const double crRightHandSideBoundedVector2 =             N(2)*data.m_ext(2);
const double crRightHandSideBoundedVector3 =             N(3)*data.m_ext(3);
const double crRightHandSideBoundedVector4 =             crRightHandSideBoundedVector0 + crRightHandSideBoundedVector1 + crRightHandSideBoundedVector2 + crRightHandSideBoundedVector3;
const double crRightHandSideBoundedVector5 =             N(0)*crRightHandSideBoundedVector4;
const double crRightHandSideBoundedVector6 =             DN_DX(0,1)*data.U(0,2);
const double crRightHandSideBoundedVector7 =             DN_DX(1,1)*data.U(1,2);
const double crRightHandSideBoundedVector8 =             DN_DX(2,1)*data.U(2,2);
const double crRightHandSideBoundedVector9 =             DN_DX(3,1)*data.U(3,2);
const double crRightHandSideBoundedVector10 =             crRightHandSideBoundedVector6 + crRightHandSideBoundedVector7 + crRightHandSideBoundedVector8 + crRightHandSideBoundedVector9;
const double crRightHandSideBoundedVector11 =             DN_DX(0,0)*data.U(0,1);
const double crRightHandSideBoundedVector12 =             DN_DX(1,0)*data.U(1,1);
const double crRightHandSideBoundedVector13 =             DN_DX(2,0)*data.U(2,1);
const double crRightHandSideBoundedVector14 =             DN_DX(3,0)*data.U(3,1);
const double crRightHandSideBoundedVector15 =             crRightHandSideBoundedVector11 + crRightHandSideBoundedVector12 + crRightHandSideBoundedVector13 + crRightHandSideBoundedVector14;
const double crRightHandSideBoundedVector16 =             crRightHandSideBoundedVector10 + crRightHandSideBoundedVector15;
const double crRightHandSideBoundedVector17 =             DN_DX(0,0)*data.U(0,0) + DN_DX(1,0)*data.U(1,0) + DN_DX(2,0)*data.U(2,0) + DN_DX(3,0)*data.U(3,0);
const double crRightHandSideBoundedVector18 =             N(0)*data.alpha_sc_nodes(0) + N(1)*data.alpha_sc_nodes(1) + N(2)*data.alpha_sc_nodes(2) + N(3)*data.alpha_sc_nodes(3);
const double crRightHandSideBoundedVector19 =             crRightHandSideBoundedVector17*crRightHandSideBoundedVector18;
const double crRightHandSideBoundedVector20 =             DN_DX(0,1)*data.U(0,0) + DN_DX(1,1)*data.U(1,0) + DN_DX(2,1)*data.U(2,0) + DN_DX(3,1)*data.U(3,0);
const double crRightHandSideBoundedVector21 =             crRightHandSideBoundedVector18*crRightHandSideBoundedVector20;
const double crRightHandSideBoundedVector22 =             N(0)*data.U(0,0) + N(1)*data.U(1,0) + N(2)*data.U(2,0) + N(3)*data.U(3,0);
const double crRightHandSideBoundedVector23 =             1.0/crRightHandSideBoundedVector22;
const double crRightHandSideBoundedVector24 =             fabs(crRightHandSideBoundedVector22);
const double crRightHandSideBoundedVector25 =             N(0)*data.U(0,1) + N(1)*data.U(1,1) + N(2)*data.U(2,1) + N(3)*data.U(3,1);
const double crRightHandSideBoundedVector26 =             pow(crRightHandSideBoundedVector25, 2);
const double crRightHandSideBoundedVector27 =             N(0)*data.U(0,2) + N(1)*data.U(1,2) + N(2)*data.U(2,2) + N(3)*data.U(3,2);
const double crRightHandSideBoundedVector28 =             pow(crRightHandSideBoundedVector27, 2);
const double crRightHandSideBoundedVector29 =             crRightHandSideBoundedVector26 + crRightHandSideBoundedVector28;
const double crRightHandSideBoundedVector30 =             0.5*crRightHandSideBoundedVector26;
const double crRightHandSideBoundedVector31 =             0.5*crRightHandSideBoundedVector28;
const double crRightHandSideBoundedVector32 =             N(0)*data.U(0,3);
const double crRightHandSideBoundedVector33 =             N(1)*data.U(1,3);
const double crRightHandSideBoundedVector34 =             N(2)*data.U(2,3);
const double crRightHandSideBoundedVector35 =             N(3)*data.U(3,3);
const double crRightHandSideBoundedVector36 =             -crRightHandSideBoundedVector32 - crRightHandSideBoundedVector33 - crRightHandSideBoundedVector34 - crRightHandSideBoundedVector35;
const double crRightHandSideBoundedVector37 =             crRightHandSideBoundedVector23*crRightHandSideBoundedVector30 + crRightHandSideBoundedVector23*crRightHandSideBoundedVector31 + crRightHandSideBoundedVector36;
const double crRightHandSideBoundedVector38 =             data.gamma - 1;
const double crRightHandSideBoundedVector39 =             crRightHandSideBoundedVector23*crRightHandSideBoundedVector38;
const double crRightHandSideBoundedVector40 =             crRightHandSideBoundedVector37*crRightHandSideBoundedVector39;
const double crRightHandSideBoundedVector41 =             1.0/data.gamma;
const double crRightHandSideBoundedVector42 =             N(0)*data.r_ext(0) + N(1)*data.r_ext(1) + N(2)*data.r_ext(2) + N(3)*data.r_ext(3);
const double crRightHandSideBoundedVector43 =             pow(crRightHandSideBoundedVector42, 2);
const double crRightHandSideBoundedVector44 =             N(0)*data.f_ext(0,0) + N(1)*data.f_ext(1,0) + N(2)*data.f_ext(2,0) + N(3)*data.f_ext(3,0);
const double crRightHandSideBoundedVector45 =             N(0)*data.f_ext(0,1) + N(1)*data.f_ext(1,1) + N(2)*data.f_ext(2,1) + N(3)*data.f_ext(3,1);
const double crRightHandSideBoundedVector46 =             crRightHandSideBoundedVector40*data.gamma*(pow(crRightHandSideBoundedVector44, 2) + pow(crRightHandSideBoundedVector45, 2));
const double crRightHandSideBoundedVector47 =             0.70710678118654757*crRightHandSideBoundedVector24*crRightHandSideBoundedVector41*stab_c3*sqrt((crRightHandSideBoundedVector43 - 2.0*crRightHandSideBoundedVector46 + 2.0*sqrt(0.25*crRightHandSideBoundedVector43 - crRightHandSideBoundedVector46)*fabs(crRightHandSideBoundedVector42))/(pow(crRightHandSideBoundedVector37, 2)*pow(crRightHandSideBoundedVector38, 2))) + stab_c2*(sqrt(data.gamma)*sqrt(-crRightHandSideBoundedVector40) + sqrt(crRightHandSideBoundedVector29)/crRightHandSideBoundedVector24)/data.h;
const double crRightHandSideBoundedVector48 =             1.0*(N(0)*data.dUdt(0,0) + N(1)*data.dUdt(1,0) + N(2)*data.dUdt(2,0) + N(3)*data.dUdt(3,0) - crRightHandSideBoundedVector0 - crRightHandSideBoundedVector1 + crRightHandSideBoundedVector16 - crRightHandSideBoundedVector2 - crRightHandSideBoundedVector3)/crRightHandSideBoundedVector47;
const double crRightHandSideBoundedVector49 =             crRightHandSideBoundedVector23*crRightHandSideBoundedVector48;
const double crRightHandSideBoundedVector50 =             crRightHandSideBoundedVector22*crRightHandSideBoundedVector44;
const double crRightHandSideBoundedVector51 =             crRightHandSideBoundedVector15*crRightHandSideBoundedVector25;
const double crRightHandSideBoundedVector52 =             1.0*data.gamma;
const double crRightHandSideBoundedVector53 =             crRightHandSideBoundedVector23*(3.0 - crRightHandSideBoundedVector52);
const double crRightHandSideBoundedVector54 =             DN_DX(0,0)*data.U(0,2);
const double crRightHandSideBoundedVector55 =             DN_DX(1,0)*data.U(1,2);
const double crRightHandSideBoundedVector56 =             DN_DX(2,0)*data.U(2,2);
const double crRightHandSideBoundedVector57 =             DN_DX(3,0)*data.U(3,2);
const double crRightHandSideBoundedVector58 =             crRightHandSideBoundedVector54 + crRightHandSideBoundedVector55 + crRightHandSideBoundedVector56 + crRightHandSideBoundedVector57;
const double crRightHandSideBoundedVector59 =             crRightHandSideBoundedVector27*crRightHandSideBoundedVector58;
const double crRightHandSideBoundedVector60 =             1.0*crRightHandSideBoundedVector39;
const double crRightHandSideBoundedVector61 =             -crRightHandSideBoundedVector59*crRightHandSideBoundedVector60;
const double crRightHandSideBoundedVector62 =             DN_DX(0,0)*data.U(0,3) + DN_DX(1,0)*data.U(1,3) + DN_DX(2,0)*data.U(2,3) + DN_DX(3,0)*data.U(3,3);
const double crRightHandSideBoundedVector63 =             DN_DX(0,1)*data.U(0,1);
const double crRightHandSideBoundedVector64 =             DN_DX(1,1)*data.U(1,1);
const double crRightHandSideBoundedVector65 =             DN_DX(2,1)*data.U(2,1);
const double crRightHandSideBoundedVector66 =             DN_DX(3,1)*data.U(3,1);
const double crRightHandSideBoundedVector67 =             crRightHandSideBoundedVector63 + crRightHandSideBoundedVector64 + crRightHandSideBoundedVector65 + crRightHandSideBoundedVector66;
const double crRightHandSideBoundedVector68 =             crRightHandSideBoundedVector27*crRightHandSideBoundedVector67;
const double crRightHandSideBoundedVector69 =             crRightHandSideBoundedVector10*crRightHandSideBoundedVector25;
const double crRightHandSideBoundedVector70 =             pow(crRightHandSideBoundedVector22, -2);
const double crRightHandSideBoundedVector71 =             0.5*crRightHandSideBoundedVector29;
const double crRightHandSideBoundedVector72 =             crRightHandSideBoundedVector38*crRightHandSideBoundedVector71;
const double crRightHandSideBoundedVector73 =             -crRightHandSideBoundedVector26 + crRightHandSideBoundedVector72;
const double crRightHandSideBoundedVector74 =             crRightHandSideBoundedVector70*crRightHandSideBoundedVector73;
const double crRightHandSideBoundedVector75 =             crRightHandSideBoundedVector20*crRightHandSideBoundedVector27;
const double crRightHandSideBoundedVector76 =             crRightHandSideBoundedVector70*crRightHandSideBoundedVector75;
const double crRightHandSideBoundedVector77 =             crRightHandSideBoundedVector17*crRightHandSideBoundedVector74 + crRightHandSideBoundedVector23*crRightHandSideBoundedVector68 + crRightHandSideBoundedVector23*crRightHandSideBoundedVector69 - crRightHandSideBoundedVector25*crRightHandSideBoundedVector76 + crRightHandSideBoundedVector38*crRightHandSideBoundedVector62 + crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector78 =             N(0)*data.mu_sc_nodes(0);
const double crRightHandSideBoundedVector79 =             N(1)*data.mu_sc_nodes(1);
const double crRightHandSideBoundedVector80 =             N(2)*data.mu_sc_nodes(2);
const double crRightHandSideBoundedVector81 =             N(3)*data.mu_sc_nodes(3);
const double crRightHandSideBoundedVector82 =             crRightHandSideBoundedVector78 + crRightHandSideBoundedVector79 + crRightHandSideBoundedVector80 + crRightHandSideBoundedVector81 + data.mu;
const double crRightHandSideBoundedVector83 =             crRightHandSideBoundedVector23*stab_c1/pow(data.h, 2);
const double crRightHandSideBoundedVector84 =             1.0/(crRightHandSideBoundedVector47 + 1.3333333333333333*crRightHandSideBoundedVector82*crRightHandSideBoundedVector83);
const double crRightHandSideBoundedVector85 =             crRightHandSideBoundedVector84*(N(0)*data.dUdt(0,1) + N(1)*data.dUdt(1,1) + N(2)*data.dUdt(2,1) + N(3)*data.dUdt(3,1) - crRightHandSideBoundedVector50 + crRightHandSideBoundedVector51*crRightHandSideBoundedVector53 + crRightHandSideBoundedVector77);
const double crRightHandSideBoundedVector86 =             crRightHandSideBoundedVector22*crRightHandSideBoundedVector45;
const double crRightHandSideBoundedVector87 =             crRightHandSideBoundedVector10*crRightHandSideBoundedVector27;
const double crRightHandSideBoundedVector88 =             crRightHandSideBoundedVector25*crRightHandSideBoundedVector67;
const double crRightHandSideBoundedVector89 =             -crRightHandSideBoundedVector60*crRightHandSideBoundedVector88;
const double crRightHandSideBoundedVector90 =             DN_DX(0,1)*data.U(0,3) + DN_DX(1,1)*data.U(1,3) + DN_DX(2,1)*data.U(2,3) + DN_DX(3,1)*data.U(3,3);
const double crRightHandSideBoundedVector91 =             crRightHandSideBoundedVector15*crRightHandSideBoundedVector27;
const double crRightHandSideBoundedVector92 =             crRightHandSideBoundedVector25*crRightHandSideBoundedVector58;
const double crRightHandSideBoundedVector93 =             -crRightHandSideBoundedVector28 + crRightHandSideBoundedVector72;
const double crRightHandSideBoundedVector94 =             crRightHandSideBoundedVector70*crRightHandSideBoundedVector93;
const double crRightHandSideBoundedVector95 =             crRightHandSideBoundedVector17*crRightHandSideBoundedVector25;
const double crRightHandSideBoundedVector96 =             crRightHandSideBoundedVector27*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector97 =             crRightHandSideBoundedVector20*crRightHandSideBoundedVector94 + crRightHandSideBoundedVector23*crRightHandSideBoundedVector91 + crRightHandSideBoundedVector23*crRightHandSideBoundedVector92 + crRightHandSideBoundedVector38*crRightHandSideBoundedVector90 + crRightHandSideBoundedVector89 - crRightHandSideBoundedVector95*crRightHandSideBoundedVector96;
const double crRightHandSideBoundedVector98 =             crRightHandSideBoundedVector84*(N(0)*data.dUdt(0,2) + N(1)*data.dUdt(1,2) + N(2)*data.dUdt(2,2) + N(3)*data.dUdt(3,2) + crRightHandSideBoundedVector53*crRightHandSideBoundedVector87 - crRightHandSideBoundedVector86 + crRightHandSideBoundedVector97);
const double crRightHandSideBoundedVector99 =             DN_DX(0,1)*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector100 =             crRightHandSideBoundedVector17*crRightHandSideBoundedVector27;
const double crRightHandSideBoundedVector101 =             crRightHandSideBoundedVector20*crRightHandSideBoundedVector25;
const double crRightHandSideBoundedVector102 =             crRightHandSideBoundedVector82*(-crRightHandSideBoundedVector100 - crRightHandSideBoundedVector101 + crRightHandSideBoundedVector22*crRightHandSideBoundedVector58 + crRightHandSideBoundedVector22*crRightHandSideBoundedVector67);
const double crRightHandSideBoundedVector103 =             crRightHandSideBoundedVector15*crRightHandSideBoundedVector22 - crRightHandSideBoundedVector95;
const double crRightHandSideBoundedVector104 =             2*crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector105 =             crRightHandSideBoundedVector10*crRightHandSideBoundedVector22 - crRightHandSideBoundedVector75;
const double crRightHandSideBoundedVector106 =             (crRightHandSideBoundedVector103 + crRightHandSideBoundedVector105)*(N(0)*data.beta_sc_nodes(0) + N(1)*data.beta_sc_nodes(1) + N(2)*data.beta_sc_nodes(2) + N(3)*data.beta_sc_nodes(3) - 0.66666666666666663*crRightHandSideBoundedVector78 - 0.66666666666666663*crRightHandSideBoundedVector79 - 0.66666666666666663*crRightHandSideBoundedVector80 - 0.66666666666666663*crRightHandSideBoundedVector81 - 0.66666666666666663*data.mu);
const double crRightHandSideBoundedVector107 =             crRightHandSideBoundedVector103*crRightHandSideBoundedVector104 + crRightHandSideBoundedVector106;
const double crRightHandSideBoundedVector108 =             DN_DX(0,0)*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector109 =             crRightHandSideBoundedVector52 - 3.0;
const double crRightHandSideBoundedVector110 =             crRightHandSideBoundedVector109*crRightHandSideBoundedVector51;
const double crRightHandSideBoundedVector111 =             -crRightHandSideBoundedVector110*crRightHandSideBoundedVector23 + crRightHandSideBoundedVector77;
const double crRightHandSideBoundedVector112 =             N(0)*crRightHandSideBoundedVector44;
const double crRightHandSideBoundedVector113 =             DN_DX(0,1)*crRightHandSideBoundedVector25;
const double crRightHandSideBoundedVector114 =             crRightHandSideBoundedVector113*crRightHandSideBoundedVector96;
const double crRightHandSideBoundedVector115 =             crRightHandSideBoundedVector23*crRightHandSideBoundedVector75;
const double crRightHandSideBoundedVector116 =             2*crRightHandSideBoundedVector115;
const double crRightHandSideBoundedVector117 =             -crRightHandSideBoundedVector116*crRightHandSideBoundedVector25 + crRightHandSideBoundedVector68 + crRightHandSideBoundedVector69;
const double crRightHandSideBoundedVector118 =             N(0)*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector119 =             1.0*crRightHandSideBoundedVector38;
const double crRightHandSideBoundedVector120 =             crRightHandSideBoundedVector17*crRightHandSideBoundedVector23;
const double crRightHandSideBoundedVector121 =             crRightHandSideBoundedVector110 + crRightHandSideBoundedVector119*crRightHandSideBoundedVector59 - 2*crRightHandSideBoundedVector120*crRightHandSideBoundedVector73;
const double crRightHandSideBoundedVector122 =             DN_DX(0,0)*crRightHandSideBoundedVector27;
const double crRightHandSideBoundedVector123 =             crRightHandSideBoundedVector101*crRightHandSideBoundedVector23 - crRightHandSideBoundedVector63 - crRightHandSideBoundedVector64 - crRightHandSideBoundedVector65 - crRightHandSideBoundedVector66;
const double crRightHandSideBoundedVector124 =             N(0)*crRightHandSideBoundedVector123;
const double crRightHandSideBoundedVector125 =             crRightHandSideBoundedVector52 - 1.0;
const double crRightHandSideBoundedVector126 =             crRightHandSideBoundedVector100*crRightHandSideBoundedVector23 - crRightHandSideBoundedVector54 - crRightHandSideBoundedVector55 - crRightHandSideBoundedVector56 - crRightHandSideBoundedVector57;
const double crRightHandSideBoundedVector127 =             N(0)*crRightHandSideBoundedVector126;
const double crRightHandSideBoundedVector128 =             crRightHandSideBoundedVector23*crRightHandSideBoundedVector98;
const double crRightHandSideBoundedVector129 =             crRightHandSideBoundedVector23*crRightHandSideBoundedVector95;
const double crRightHandSideBoundedVector130 =             -crRightHandSideBoundedVector11 - crRightHandSideBoundedVector12 + crRightHandSideBoundedVector129 - crRightHandSideBoundedVector13 - crRightHandSideBoundedVector14;
const double crRightHandSideBoundedVector131 =             N(0)*crRightHandSideBoundedVector130;
const double crRightHandSideBoundedVector132 =             DN_DX(0,0)*crRightHandSideBoundedVector25;
const double crRightHandSideBoundedVector133 =             DN_DX(0,1)*crRightHandSideBoundedVector27;
const double crRightHandSideBoundedVector134 =             crRightHandSideBoundedVector115 - crRightHandSideBoundedVector6 - crRightHandSideBoundedVector7 - crRightHandSideBoundedVector8 - crRightHandSideBoundedVector9;
const double crRightHandSideBoundedVector135 =             N(0)*crRightHandSideBoundedVector134;
const double crRightHandSideBoundedVector136 =             crRightHandSideBoundedVector133 - crRightHandSideBoundedVector135;
const double crRightHandSideBoundedVector137 =             crRightHandSideBoundedVector23*crRightHandSideBoundedVector85;
const double crRightHandSideBoundedVector138 =             (N(0)*data.lamb_sc_nodes(0) + N(1)*data.lamb_sc_nodes(1) + N(2)*data.lamb_sc_nodes(2) + N(3)*data.lamb_sc_nodes(3) + data.lambda)/data.c_v;
const double crRightHandSideBoundedVector139 =             crRightHandSideBoundedVector22*crRightHandSideBoundedVector42;
const double crRightHandSideBoundedVector140 =             crRightHandSideBoundedVector25*crRightHandSideBoundedVector44;
const double crRightHandSideBoundedVector141 =             crRightHandSideBoundedVector27*crRightHandSideBoundedVector45;
const double crRightHandSideBoundedVector142 =             crRightHandSideBoundedVector62*data.gamma;
const double crRightHandSideBoundedVector143 =             crRightHandSideBoundedVector23*crRightHandSideBoundedVector25;
const double crRightHandSideBoundedVector144 =             crRightHandSideBoundedVector90*data.gamma;
const double crRightHandSideBoundedVector145 =             crRightHandSideBoundedVector23*crRightHandSideBoundedVector27;
const double crRightHandSideBoundedVector146 =             1.0*crRightHandSideBoundedVector92;
const double crRightHandSideBoundedVector147 =             1.0*crRightHandSideBoundedVector68;
const double crRightHandSideBoundedVector148 =             crRightHandSideBoundedVector147*crRightHandSideBoundedVector25;
const double crRightHandSideBoundedVector149 =             crRightHandSideBoundedVector26*crRightHandSideBoundedVector60;
const double crRightHandSideBoundedVector150 =             crRightHandSideBoundedVector32 + crRightHandSideBoundedVector33 + crRightHandSideBoundedVector34 + crRightHandSideBoundedVector35;
const double crRightHandSideBoundedVector151 =             crRightHandSideBoundedVector38*(crRightHandSideBoundedVector150 - crRightHandSideBoundedVector23*(crRightHandSideBoundedVector30 + crRightHandSideBoundedVector31));
const double crRightHandSideBoundedVector152 =             crRightHandSideBoundedVector150 + crRightHandSideBoundedVector151;
const double crRightHandSideBoundedVector153 =             crRightHandSideBoundedVector28*crRightHandSideBoundedVector60;
const double crRightHandSideBoundedVector154 =             -crRightHandSideBoundedVector151 + crRightHandSideBoundedVector36;
const double crRightHandSideBoundedVector155 =             crRightHandSideBoundedVector154 + crRightHandSideBoundedVector39*crRightHandSideBoundedVector71;
const double crRightHandSideBoundedVector156 =             crRightHandSideBoundedVector155*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector157 =             (N(0)*data.dUdt(0,3) + N(1)*data.dUdt(1,3) + N(2)*data.dUdt(2,3) + N(3)*data.dUdt(3,3) + crRightHandSideBoundedVector10*crRightHandSideBoundedVector23*(crRightHandSideBoundedVector152 - crRightHandSideBoundedVector153) - crRightHandSideBoundedVector139 - crRightHandSideBoundedVector140 - crRightHandSideBoundedVector141 + crRightHandSideBoundedVector142*crRightHandSideBoundedVector143 + crRightHandSideBoundedVector144*crRightHandSideBoundedVector145 - crRightHandSideBoundedVector146*crRightHandSideBoundedVector38*crRightHandSideBoundedVector96 - crRightHandSideBoundedVector148*crRightHandSideBoundedVector38*crRightHandSideBoundedVector70 + crRightHandSideBoundedVector15*crRightHandSideBoundedVector23*(-crRightHandSideBoundedVector149 + crRightHandSideBoundedVector152) + crRightHandSideBoundedVector155*crRightHandSideBoundedVector76 + crRightHandSideBoundedVector156*crRightHandSideBoundedVector95)/(crRightHandSideBoundedVector138*crRightHandSideBoundedVector41*crRightHandSideBoundedVector83 + crRightHandSideBoundedVector47);
const double crRightHandSideBoundedVector158 =             crRightHandSideBoundedVector119*crRightHandSideBoundedVector157;
const double crRightHandSideBoundedVector159 =             crRightHandSideBoundedVector104*crRightHandSideBoundedVector105 + crRightHandSideBoundedVector106;
const double crRightHandSideBoundedVector160 =             crRightHandSideBoundedVector109*crRightHandSideBoundedVector87;
const double crRightHandSideBoundedVector161 =             -crRightHandSideBoundedVector160*crRightHandSideBoundedVector23 + crRightHandSideBoundedVector97;
const double crRightHandSideBoundedVector162 =             N(0)*crRightHandSideBoundedVector45;
const double crRightHandSideBoundedVector163 =             crRightHandSideBoundedVector132*crRightHandSideBoundedVector96;
const double crRightHandSideBoundedVector164 =             2*crRightHandSideBoundedVector129;
const double crRightHandSideBoundedVector165 =             -crRightHandSideBoundedVector164*crRightHandSideBoundedVector27 + crRightHandSideBoundedVector91 + crRightHandSideBoundedVector92;
const double crRightHandSideBoundedVector166 =             crRightHandSideBoundedVector20*crRightHandSideBoundedVector23;
const double crRightHandSideBoundedVector167 =             crRightHandSideBoundedVector119*crRightHandSideBoundedVector88 + crRightHandSideBoundedVector160 - 2*crRightHandSideBoundedVector166*crRightHandSideBoundedVector93;
const double crRightHandSideBoundedVector168 =             -crRightHandSideBoundedVector131 + crRightHandSideBoundedVector132;
const double crRightHandSideBoundedVector169 =             crRightHandSideBoundedVector139 + crRightHandSideBoundedVector140 + crRightHandSideBoundedVector141;
const double crRightHandSideBoundedVector170 =             crRightHandSideBoundedVector102*crRightHandSideBoundedVector145 + crRightHandSideBoundedVector107*crRightHandSideBoundedVector143 + crRightHandSideBoundedVector138*(crRightHandSideBoundedVector120*crRightHandSideBoundedVector26 + crRightHandSideBoundedVector120*crRightHandSideBoundedVector28 - crRightHandSideBoundedVector150*crRightHandSideBoundedVector17 + crRightHandSideBoundedVector22*crRightHandSideBoundedVector62 - crRightHandSideBoundedVector51 - crRightHandSideBoundedVector59);
const double crRightHandSideBoundedVector171 =             crRightHandSideBoundedVector102*crRightHandSideBoundedVector143 + crRightHandSideBoundedVector138*(-crRightHandSideBoundedVector150*crRightHandSideBoundedVector20 + crRightHandSideBoundedVector166*crRightHandSideBoundedVector26 + crRightHandSideBoundedVector166*crRightHandSideBoundedVector28 + crRightHandSideBoundedVector22*crRightHandSideBoundedVector90 - crRightHandSideBoundedVector87 - crRightHandSideBoundedVector88) + crRightHandSideBoundedVector145*crRightHandSideBoundedVector159;
const double crRightHandSideBoundedVector172 =             crRightHandSideBoundedVector142*crRightHandSideBoundedVector25;
const double crRightHandSideBoundedVector173 =             crRightHandSideBoundedVector144*crRightHandSideBoundedVector27;
const double crRightHandSideBoundedVector174 =             crRightHandSideBoundedVector149 + crRightHandSideBoundedVector154;
const double crRightHandSideBoundedVector175 =             crRightHandSideBoundedVector153 + crRightHandSideBoundedVector154;
const double crRightHandSideBoundedVector176 =             -crRightHandSideBoundedVector10*crRightHandSideBoundedVector175 + crRightHandSideBoundedVector115*crRightHandSideBoundedVector155 + crRightHandSideBoundedVector129*crRightHandSideBoundedVector155 - crRightHandSideBoundedVector146*crRightHandSideBoundedVector27*crRightHandSideBoundedVector39 - crRightHandSideBoundedVector148*crRightHandSideBoundedVector39 - crRightHandSideBoundedVector15*crRightHandSideBoundedVector174 + crRightHandSideBoundedVector172 + crRightHandSideBoundedVector173;
const double crRightHandSideBoundedVector177 =             N(0)*crRightHandSideBoundedVector23;
const double crRightHandSideBoundedVector178 =             -2.0*crRightHandSideBoundedVector115*crRightHandSideBoundedVector25 + crRightHandSideBoundedVector147 + 1.0*crRightHandSideBoundedVector69;
const double crRightHandSideBoundedVector179 =             crRightHandSideBoundedVector118*crRightHandSideBoundedVector38;
const double crRightHandSideBoundedVector180 =             crRightHandSideBoundedVector174*crRightHandSideBoundedVector23;
const double crRightHandSideBoundedVector181 =             3.0*crRightHandSideBoundedVector39;
const double crRightHandSideBoundedVector182 =             2.0*crRightHandSideBoundedVector39;
const double crRightHandSideBoundedVector183 =             crRightHandSideBoundedVector155 + crRightHandSideBoundedVector182*crRightHandSideBoundedVector26;
const double crRightHandSideBoundedVector184 =             crRightHandSideBoundedVector120*crRightHandSideBoundedVector183 + crRightHandSideBoundedVector142 - crRightHandSideBoundedVector181*crRightHandSideBoundedVector51 + crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector185 =             -2.0*crRightHandSideBoundedVector129*crRightHandSideBoundedVector27 + crRightHandSideBoundedVector146 + 1.0*crRightHandSideBoundedVector91;
const double crRightHandSideBoundedVector186 =             crRightHandSideBoundedVector175*crRightHandSideBoundedVector23;
const double crRightHandSideBoundedVector187 =             crRightHandSideBoundedVector155 + crRightHandSideBoundedVector182*crRightHandSideBoundedVector28;
const double crRightHandSideBoundedVector188 =             crRightHandSideBoundedVector144 + crRightHandSideBoundedVector166*crRightHandSideBoundedVector187 - crRightHandSideBoundedVector181*crRightHandSideBoundedVector87 + crRightHandSideBoundedVector89;
const double crRightHandSideBoundedVector189 =             crRightHandSideBoundedVector157*crRightHandSideBoundedVector23*crRightHandSideBoundedVector52;
const double crRightHandSideBoundedVector190 =             crRightHandSideBoundedVector154 + crRightHandSideBoundedVector29*crRightHandSideBoundedVector39;
const double crRightHandSideBoundedVector191 =             crRightHandSideBoundedVector15*crRightHandSideBoundedVector183 - crRightHandSideBoundedVector164*crRightHandSideBoundedVector190 - crRightHandSideBoundedVector172 + crRightHandSideBoundedVector182*crRightHandSideBoundedVector27*crRightHandSideBoundedVector92;
const double crRightHandSideBoundedVector192 =             crRightHandSideBoundedVector10*crRightHandSideBoundedVector187 - crRightHandSideBoundedVector116*crRightHandSideBoundedVector190 - crRightHandSideBoundedVector173 + crRightHandSideBoundedVector182*crRightHandSideBoundedVector25*crRightHandSideBoundedVector68;
const double crRightHandSideBoundedVector193 =             N(1)*crRightHandSideBoundedVector4;
const double crRightHandSideBoundedVector194 =             DN_DX(1,1)*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector195 =             DN_DX(1,0)*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector196 =             N(1)*crRightHandSideBoundedVector44;
const double crRightHandSideBoundedVector197 =             DN_DX(1,1)*crRightHandSideBoundedVector25;
const double crRightHandSideBoundedVector198 =             crRightHandSideBoundedVector197*crRightHandSideBoundedVector96;
const double crRightHandSideBoundedVector199 =             N(1)*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector200 =             DN_DX(1,0)*crRightHandSideBoundedVector27;
const double crRightHandSideBoundedVector201 =             N(1)*crRightHandSideBoundedVector123;
const double crRightHandSideBoundedVector202 =             N(1)*crRightHandSideBoundedVector126;
const double crRightHandSideBoundedVector203 =             N(1)*crRightHandSideBoundedVector130;
const double crRightHandSideBoundedVector204 =             DN_DX(1,0)*crRightHandSideBoundedVector25;
const double crRightHandSideBoundedVector205 =             DN_DX(1,1)*crRightHandSideBoundedVector27;
const double crRightHandSideBoundedVector206 =             N(1)*crRightHandSideBoundedVector134;
const double crRightHandSideBoundedVector207 =             crRightHandSideBoundedVector205 - crRightHandSideBoundedVector206;
const double crRightHandSideBoundedVector208 =             N(1)*crRightHandSideBoundedVector45;
const double crRightHandSideBoundedVector209 =             crRightHandSideBoundedVector204*crRightHandSideBoundedVector96;
const double crRightHandSideBoundedVector210 =             -crRightHandSideBoundedVector203 + crRightHandSideBoundedVector204;
const double crRightHandSideBoundedVector211 =             N(1)*crRightHandSideBoundedVector23;
const double crRightHandSideBoundedVector212 =             crRightHandSideBoundedVector199*crRightHandSideBoundedVector38;
const double crRightHandSideBoundedVector213 =             N(2)*crRightHandSideBoundedVector4;
const double crRightHandSideBoundedVector214 =             DN_DX(2,1)*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector215 =             DN_DX(2,0)*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector216 =             N(2)*crRightHandSideBoundedVector44;
const double crRightHandSideBoundedVector217 =             DN_DX(2,1)*crRightHandSideBoundedVector25;
const double crRightHandSideBoundedVector218 =             crRightHandSideBoundedVector217*crRightHandSideBoundedVector96;
const double crRightHandSideBoundedVector219 =             N(2)*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector220 =             DN_DX(2,0)*crRightHandSideBoundedVector27;
const double crRightHandSideBoundedVector221 =             N(2)*crRightHandSideBoundedVector123;
const double crRightHandSideBoundedVector222 =             N(2)*crRightHandSideBoundedVector126;
const double crRightHandSideBoundedVector223 =             N(2)*crRightHandSideBoundedVector130;
const double crRightHandSideBoundedVector224 =             DN_DX(2,0)*crRightHandSideBoundedVector25;
const double crRightHandSideBoundedVector225 =             DN_DX(2,1)*crRightHandSideBoundedVector27;
const double crRightHandSideBoundedVector226 =             N(2)*crRightHandSideBoundedVector134;
const double crRightHandSideBoundedVector227 =             crRightHandSideBoundedVector225 - crRightHandSideBoundedVector226;
const double crRightHandSideBoundedVector228 =             N(2)*crRightHandSideBoundedVector45;
const double crRightHandSideBoundedVector229 =             crRightHandSideBoundedVector224*crRightHandSideBoundedVector96;
const double crRightHandSideBoundedVector230 =             -crRightHandSideBoundedVector223 + crRightHandSideBoundedVector224;
const double crRightHandSideBoundedVector231 =             N(2)*crRightHandSideBoundedVector23;
const double crRightHandSideBoundedVector232 =             crRightHandSideBoundedVector219*crRightHandSideBoundedVector38;
const double crRightHandSideBoundedVector233 =             N(3)*crRightHandSideBoundedVector4;
const double crRightHandSideBoundedVector234 =             DN_DX(3,1)*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector235 =             DN_DX(3,0)*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector236 =             N(3)*crRightHandSideBoundedVector44;
const double crRightHandSideBoundedVector237 =             DN_DX(3,1)*crRightHandSideBoundedVector25;
const double crRightHandSideBoundedVector238 =             crRightHandSideBoundedVector237*crRightHandSideBoundedVector96;
const double crRightHandSideBoundedVector239 =             N(3)*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector240 =             DN_DX(3,0)*crRightHandSideBoundedVector27;
const double crRightHandSideBoundedVector241 =             N(3)*crRightHandSideBoundedVector123;
const double crRightHandSideBoundedVector242 =             N(3)*crRightHandSideBoundedVector126;
const double crRightHandSideBoundedVector243 =             N(3)*crRightHandSideBoundedVector130;
const double crRightHandSideBoundedVector244 =             DN_DX(3,0)*crRightHandSideBoundedVector25;
const double crRightHandSideBoundedVector245 =             DN_DX(3,1)*crRightHandSideBoundedVector27;
const double crRightHandSideBoundedVector246 =             N(3)*crRightHandSideBoundedVector134;
const double crRightHandSideBoundedVector247 =             crRightHandSideBoundedVector245 - crRightHandSideBoundedVector246;
const double crRightHandSideBoundedVector248 =             N(3)*crRightHandSideBoundedVector45;
const double crRightHandSideBoundedVector249 =             crRightHandSideBoundedVector244*crRightHandSideBoundedVector96;
const double crRightHandSideBoundedVector250 =             -crRightHandSideBoundedVector243 + crRightHandSideBoundedVector244;
const double crRightHandSideBoundedVector251 =             N(3)*crRightHandSideBoundedVector23;
const double crRightHandSideBoundedVector252 =             crRightHandSideBoundedVector239*crRightHandSideBoundedVector38;
            rRightHandSideBoundedVector[0] += DN_DX(0,0)*crRightHandSideBoundedVector19 - DN_DX(0,0)*crRightHandSideBoundedVector85 + DN_DX(0,1)*crRightHandSideBoundedVector21 - DN_DX(0,1)*crRightHandSideBoundedVector98 - N(0)*crRightHandSideBoundedVector16 - crRightHandSideBoundedVector49*crRightHandSideBoundedVector5 + crRightHandSideBoundedVector5;
            rRightHandSideBoundedVector[1] += -DN_DX(0,0)*crRightHandSideBoundedVector158 - N(0)*crRightHandSideBoundedVector111 + N(0)*crRightHandSideBoundedVector50 - crRightHandSideBoundedVector102*crRightHandSideBoundedVector99 - crRightHandSideBoundedVector107*crRightHandSideBoundedVector108 - crRightHandSideBoundedVector128*(crRightHandSideBoundedVector113 - crRightHandSideBoundedVector119*crRightHandSideBoundedVector122 - crRightHandSideBoundedVector124 + crRightHandSideBoundedVector125*crRightHandSideBoundedVector127) - crRightHandSideBoundedVector137*(crRightHandSideBoundedVector109*crRightHandSideBoundedVector131 - crRightHandSideBoundedVector109*crRightHandSideBoundedVector132 + crRightHandSideBoundedVector136) - crRightHandSideBoundedVector48*(DN_DX(0,0)*crRightHandSideBoundedVector74 + crRightHandSideBoundedVector112 - crRightHandSideBoundedVector114 - crRightHandSideBoundedVector117*crRightHandSideBoundedVector118 + crRightHandSideBoundedVector118*crRightHandSideBoundedVector121);
            rRightHandSideBoundedVector[2] += -DN_DX(0,1)*crRightHandSideBoundedVector158 - N(0)*crRightHandSideBoundedVector161 + N(0)*crRightHandSideBoundedVector86 - crRightHandSideBoundedVector102*crRightHandSideBoundedVector108 - crRightHandSideBoundedVector128*(-crRightHandSideBoundedVector109*crRightHandSideBoundedVector133 + crRightHandSideBoundedVector109*crRightHandSideBoundedVector135 + crRightHandSideBoundedVector168) - crRightHandSideBoundedVector137*(-crRightHandSideBoundedVector113*crRightHandSideBoundedVector119 + crRightHandSideBoundedVector122 + crRightHandSideBoundedVector124*crRightHandSideBoundedVector125 - crRightHandSideBoundedVector127) - crRightHandSideBoundedVector159*crRightHandSideBoundedVector99 - crRightHandSideBoundedVector48*(DN_DX(0,1)*crRightHandSideBoundedVector94 - crRightHandSideBoundedVector118*crRightHandSideBoundedVector165 + crRightHandSideBoundedVector118*crRightHandSideBoundedVector167 + crRightHandSideBoundedVector162 - crRightHandSideBoundedVector163);
            rRightHandSideBoundedVector[3] += N(0)*crRightHandSideBoundedVector169 - crRightHandSideBoundedVector108*crRightHandSideBoundedVector170 - crRightHandSideBoundedVector171*crRightHandSideBoundedVector99 - crRightHandSideBoundedVector176*crRightHandSideBoundedVector177 - crRightHandSideBoundedVector189*(crRightHandSideBoundedVector136 + crRightHandSideBoundedVector168) - crRightHandSideBoundedVector48*(N(0)*crRightHandSideBoundedVector42 + crRightHandSideBoundedVector118*crRightHandSideBoundedVector191 + crRightHandSideBoundedVector118*crRightHandSideBoundedVector192 + crRightHandSideBoundedVector132*crRightHandSideBoundedVector156 + crRightHandSideBoundedVector133*crRightHandSideBoundedVector156) - crRightHandSideBoundedVector85*(-DN_DX(0,0)*crRightHandSideBoundedVector180 + crRightHandSideBoundedVector112 - crRightHandSideBoundedVector114*crRightHandSideBoundedVector119 + crRightHandSideBoundedVector177*crRightHandSideBoundedVector184 - crRightHandSideBoundedVector178*crRightHandSideBoundedVector179) - crRightHandSideBoundedVector98*(-DN_DX(0,1)*crRightHandSideBoundedVector186 - crRightHandSideBoundedVector119*crRightHandSideBoundedVector163 + crRightHandSideBoundedVector162 + crRightHandSideBoundedVector177*crRightHandSideBoundedVector188 - crRightHandSideBoundedVector179*crRightHandSideBoundedVector185);
            rRightHandSideBoundedVector[4] += DN_DX(1,0)*crRightHandSideBoundedVector19 - DN_DX(1,0)*crRightHandSideBoundedVector85 + DN_DX(1,1)*crRightHandSideBoundedVector21 - DN_DX(1,1)*crRightHandSideBoundedVector98 - N(1)*crRightHandSideBoundedVector16 - crRightHandSideBoundedVector193*crRightHandSideBoundedVector49 + crRightHandSideBoundedVector193;
            rRightHandSideBoundedVector[5] += -DN_DX(1,0)*crRightHandSideBoundedVector158 - N(1)*crRightHandSideBoundedVector111 + N(1)*crRightHandSideBoundedVector50 - crRightHandSideBoundedVector102*crRightHandSideBoundedVector194 - crRightHandSideBoundedVector107*crRightHandSideBoundedVector195 - crRightHandSideBoundedVector128*(-crRightHandSideBoundedVector119*crRightHandSideBoundedVector200 + crRightHandSideBoundedVector125*crRightHandSideBoundedVector202 + crRightHandSideBoundedVector197 - crRightHandSideBoundedVector201) - crRightHandSideBoundedVector137*(crRightHandSideBoundedVector109*crRightHandSideBoundedVector203 - crRightHandSideBoundedVector109*crRightHandSideBoundedVector204 + crRightHandSideBoundedVector207) - crRightHandSideBoundedVector48*(DN_DX(1,0)*crRightHandSideBoundedVector74 - crRightHandSideBoundedVector117*crRightHandSideBoundedVector199 + crRightHandSideBoundedVector121*crRightHandSideBoundedVector199 + crRightHandSideBoundedVector196 - crRightHandSideBoundedVector198);
            rRightHandSideBoundedVector[6] += -DN_DX(1,1)*crRightHandSideBoundedVector158 - N(1)*crRightHandSideBoundedVector161 + N(1)*crRightHandSideBoundedVector86 - crRightHandSideBoundedVector102*crRightHandSideBoundedVector195 - crRightHandSideBoundedVector128*(-crRightHandSideBoundedVector109*crRightHandSideBoundedVector205 + crRightHandSideBoundedVector109*crRightHandSideBoundedVector206 + crRightHandSideBoundedVector210) - crRightHandSideBoundedVector137*(-crRightHandSideBoundedVector119*crRightHandSideBoundedVector197 + crRightHandSideBoundedVector125*crRightHandSideBoundedVector201 + crRightHandSideBoundedVector200 - crRightHandSideBoundedVector202) - crRightHandSideBoundedVector159*crRightHandSideBoundedVector194 - crRightHandSideBoundedVector48*(DN_DX(1,1)*crRightHandSideBoundedVector94 - crRightHandSideBoundedVector165*crRightHandSideBoundedVector199 + crRightHandSideBoundedVector167*crRightHandSideBoundedVector199 + crRightHandSideBoundedVector208 - crRightHandSideBoundedVector209);
            rRightHandSideBoundedVector[7] += N(1)*crRightHandSideBoundedVector169 - crRightHandSideBoundedVector170*crRightHandSideBoundedVector195 - crRightHandSideBoundedVector171*crRightHandSideBoundedVector194 - crRightHandSideBoundedVector176*crRightHandSideBoundedVector211 - crRightHandSideBoundedVector189*(crRightHandSideBoundedVector207 + crRightHandSideBoundedVector210) - crRightHandSideBoundedVector48*(N(1)*crRightHandSideBoundedVector42 + crRightHandSideBoundedVector156*crRightHandSideBoundedVector204 + crRightHandSideBoundedVector156*crRightHandSideBoundedVector205 + crRightHandSideBoundedVector191*crRightHandSideBoundedVector199 + crRightHandSideBoundedVector192*crRightHandSideBoundedVector199) - crRightHandSideBoundedVector85*(-DN_DX(1,0)*crRightHandSideBoundedVector180 - crRightHandSideBoundedVector119*crRightHandSideBoundedVector198 - crRightHandSideBoundedVector178*crRightHandSideBoundedVector212 + crRightHandSideBoundedVector184*crRightHandSideBoundedVector211 + crRightHandSideBoundedVector196) - crRightHandSideBoundedVector98*(-DN_DX(1,1)*crRightHandSideBoundedVector186 - crRightHandSideBoundedVector119*crRightHandSideBoundedVector209 - crRightHandSideBoundedVector185*crRightHandSideBoundedVector212 + crRightHandSideBoundedVector188*crRightHandSideBoundedVector211 + crRightHandSideBoundedVector208);
            rRightHandSideBoundedVector[8] += DN_DX(2,0)*crRightHandSideBoundedVector19 - DN_DX(2,0)*crRightHandSideBoundedVector85 + DN_DX(2,1)*crRightHandSideBoundedVector21 - DN_DX(2,1)*crRightHandSideBoundedVector98 - N(2)*crRightHandSideBoundedVector16 - crRightHandSideBoundedVector213*crRightHandSideBoundedVector49 + crRightHandSideBoundedVector213;
            rRightHandSideBoundedVector[9] += -DN_DX(2,0)*crRightHandSideBoundedVector158 - N(2)*crRightHandSideBoundedVector111 + N(2)*crRightHandSideBoundedVector50 - crRightHandSideBoundedVector102*crRightHandSideBoundedVector214 - crRightHandSideBoundedVector107*crRightHandSideBoundedVector215 - crRightHandSideBoundedVector128*(-crRightHandSideBoundedVector119*crRightHandSideBoundedVector220 + crRightHandSideBoundedVector125*crRightHandSideBoundedVector222 + crRightHandSideBoundedVector217 - crRightHandSideBoundedVector221) - crRightHandSideBoundedVector137*(crRightHandSideBoundedVector109*crRightHandSideBoundedVector223 - crRightHandSideBoundedVector109*crRightHandSideBoundedVector224 + crRightHandSideBoundedVector227) - crRightHandSideBoundedVector48*(DN_DX(2,0)*crRightHandSideBoundedVector74 - crRightHandSideBoundedVector117*crRightHandSideBoundedVector219 + crRightHandSideBoundedVector121*crRightHandSideBoundedVector219 + crRightHandSideBoundedVector216 - crRightHandSideBoundedVector218);
            rRightHandSideBoundedVector[10] += -DN_DX(2,1)*crRightHandSideBoundedVector158 - N(2)*crRightHandSideBoundedVector161 + N(2)*crRightHandSideBoundedVector86 - crRightHandSideBoundedVector102*crRightHandSideBoundedVector215 - crRightHandSideBoundedVector128*(-crRightHandSideBoundedVector109*crRightHandSideBoundedVector225 + crRightHandSideBoundedVector109*crRightHandSideBoundedVector226 + crRightHandSideBoundedVector230) - crRightHandSideBoundedVector137*(-crRightHandSideBoundedVector119*crRightHandSideBoundedVector217 + crRightHandSideBoundedVector125*crRightHandSideBoundedVector221 + crRightHandSideBoundedVector220 - crRightHandSideBoundedVector222) - crRightHandSideBoundedVector159*crRightHandSideBoundedVector214 - crRightHandSideBoundedVector48*(DN_DX(2,1)*crRightHandSideBoundedVector94 - crRightHandSideBoundedVector165*crRightHandSideBoundedVector219 + crRightHandSideBoundedVector167*crRightHandSideBoundedVector219 + crRightHandSideBoundedVector228 - crRightHandSideBoundedVector229);
            rRightHandSideBoundedVector[11] += N(2)*crRightHandSideBoundedVector169 - crRightHandSideBoundedVector170*crRightHandSideBoundedVector215 - crRightHandSideBoundedVector171*crRightHandSideBoundedVector214 - crRightHandSideBoundedVector176*crRightHandSideBoundedVector231 - crRightHandSideBoundedVector189*(crRightHandSideBoundedVector227 + crRightHandSideBoundedVector230) - crRightHandSideBoundedVector48*(N(2)*crRightHandSideBoundedVector42 + crRightHandSideBoundedVector156*crRightHandSideBoundedVector224 + crRightHandSideBoundedVector156*crRightHandSideBoundedVector225 + crRightHandSideBoundedVector191*crRightHandSideBoundedVector219 + crRightHandSideBoundedVector192*crRightHandSideBoundedVector219) - crRightHandSideBoundedVector85*(-DN_DX(2,0)*crRightHandSideBoundedVector180 - crRightHandSideBoundedVector119*crRightHandSideBoundedVector218 - crRightHandSideBoundedVector178*crRightHandSideBoundedVector232 + crRightHandSideBoundedVector184*crRightHandSideBoundedVector231 + crRightHandSideBoundedVector216) - crRightHandSideBoundedVector98*(-DN_DX(2,1)*crRightHandSideBoundedVector186 - crRightHandSideBoundedVector119*crRightHandSideBoundedVector229 - crRightHandSideBoundedVector185*crRightHandSideBoundedVector232 + crRightHandSideBoundedVector188*crRightHandSideBoundedVector231 + crRightHandSideBoundedVector228);
            rRightHandSideBoundedVector[12] += DN_DX(3,0)*crRightHandSideBoundedVector19 - DN_DX(3,0)*crRightHandSideBoundedVector85 + DN_DX(3,1)*crRightHandSideBoundedVector21 - DN_DX(3,1)*crRightHandSideBoundedVector98 - N(3)*crRightHandSideBoundedVector16 - crRightHandSideBoundedVector233*crRightHandSideBoundedVector49 + crRightHandSideBoundedVector233;
            rRightHandSideBoundedVector[13] += -DN_DX(3,0)*crRightHandSideBoundedVector158 - N(3)*crRightHandSideBoundedVector111 + N(3)*crRightHandSideBoundedVector50 - crRightHandSideBoundedVector102*crRightHandSideBoundedVector234 - crRightHandSideBoundedVector107*crRightHandSideBoundedVector235 - crRightHandSideBoundedVector128*(-crRightHandSideBoundedVector119*crRightHandSideBoundedVector240 + crRightHandSideBoundedVector125*crRightHandSideBoundedVector242 + crRightHandSideBoundedVector237 - crRightHandSideBoundedVector241) - crRightHandSideBoundedVector137*(crRightHandSideBoundedVector109*crRightHandSideBoundedVector243 - crRightHandSideBoundedVector109*crRightHandSideBoundedVector244 + crRightHandSideBoundedVector247) - crRightHandSideBoundedVector48*(DN_DX(3,0)*crRightHandSideBoundedVector74 - crRightHandSideBoundedVector117*crRightHandSideBoundedVector239 + crRightHandSideBoundedVector121*crRightHandSideBoundedVector239 + crRightHandSideBoundedVector236 - crRightHandSideBoundedVector238);
            rRightHandSideBoundedVector[14] += -DN_DX(3,1)*crRightHandSideBoundedVector158 - N(3)*crRightHandSideBoundedVector161 + N(3)*crRightHandSideBoundedVector86 - crRightHandSideBoundedVector102*crRightHandSideBoundedVector235 - crRightHandSideBoundedVector128*(-crRightHandSideBoundedVector109*crRightHandSideBoundedVector245 + crRightHandSideBoundedVector109*crRightHandSideBoundedVector246 + crRightHandSideBoundedVector250) - crRightHandSideBoundedVector137*(-crRightHandSideBoundedVector119*crRightHandSideBoundedVector237 + crRightHandSideBoundedVector125*crRightHandSideBoundedVector241 + crRightHandSideBoundedVector240 - crRightHandSideBoundedVector242) - crRightHandSideBoundedVector159*crRightHandSideBoundedVector234 - crRightHandSideBoundedVector48*(DN_DX(3,1)*crRightHandSideBoundedVector94 - crRightHandSideBoundedVector165*crRightHandSideBoundedVector239 + crRightHandSideBoundedVector167*crRightHandSideBoundedVector239 + crRightHandSideBoundedVector248 - crRightHandSideBoundedVector249);
            rRightHandSideBoundedVector[15] += N(3)*crRightHandSideBoundedVector169 - crRightHandSideBoundedVector170*crRightHandSideBoundedVector235 - crRightHandSideBoundedVector171*crRightHandSideBoundedVector234 - crRightHandSideBoundedVector176*crRightHandSideBoundedVector251 - crRightHandSideBoundedVector189*(crRightHandSideBoundedVector247 + crRightHandSideBoundedVector250) - crRightHandSideBoundedVector48*(N(3)*crRightHandSideBoundedVector42 + crRightHandSideBoundedVector156*crRightHandSideBoundedVector244 + crRightHandSideBoundedVector156*crRightHandSideBoundedVector245 + crRightHandSideBoundedVector191*crRightHandSideBoundedVector239 + crRightHandSideBoundedVector192*crRightHandSideBoundedVector239) - crRightHandSideBoundedVector85*(-DN_DX(3,0)*crRightHandSideBoundedVector180 - crRightHandSideBoundedVector119*crRightHandSideBoundedVector238 - crRightHandSideBoundedVector178*crRightHandSideBoundedVector252 + crRightHandSideBoundedVector184*crRightHandSideBoundedVector251 + crRightHandSideBoundedVector236) - crRightHandSideBoundedVector98*(-DN_DX(3,1)*crRightHandSideBoundedVector186 - crRightHandSideBoundedVector119*crRightHandSideBoundedVector249 - crRightHandSideBoundedVector185*crRightHandSideBoundedVector252 + crRightHandSideBoundedVector188*crRightHandSideBoundedVector251 + crRightHandSideBoundedVector248);

        }
    }

    KRATOS_CATCH("")
}


template <>
void CompressibleNavierStokesExplicit<2,4>::CalculateMassMatrix(
    MatrixType &rMassMatrix,
    const ProcessInfo &rCurrentProcessInfo)
{
    constexpr double dof_weight = 1.0 / DofSize;

    rMassMatrix = ZeroMatrix(DofSize, DofSize);

    for(IndexType i=0; i<NumNodes; ++i)
    {
        for(IndexType j=i; j<NumNodes; ++j)
        {
            const IndexType dof = i + j * BlockSize;

            rMassMatrix(i, dof) += dof_weight;
            rMassMatrix(dof, i) += dof_weight;
        }
    }

    // Here we assume that all the Gauss pt. have the same weight so we multiply by the volume
    rMassMatrix *= GetGeometry().Area();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class CompressibleNavierStokesExplicit<2,4>;

}
