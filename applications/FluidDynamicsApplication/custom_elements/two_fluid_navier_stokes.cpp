//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Daniel Diez
//  Co-authors:      Ruben Zorrilla
//

#include "two_fluid_navier_stokes.h"
#include "custom_utilities/two_fluid_navier_stokes_data.h"
#include "custom_utilities/two_fluid_navier_stokes_alpha_method_data.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template <class TElementData>
TwoFluidNavierStokes<TElementData>::TwoFluidNavierStokes(IndexType NewId)
    : FluidElement<TElementData>(NewId) {}

template <class TElementData>
TwoFluidNavierStokes<TElementData>::TwoFluidNavierStokes(
    IndexType NewId, const NodesArrayType &ThisNodes)
    : FluidElement<TElementData>(NewId, ThisNodes) {}

template <class TElementData>
TwoFluidNavierStokes<TElementData>::TwoFluidNavierStokes(
    IndexType NewId, GeometryType::Pointer pGeometry)
    : FluidElement<TElementData>(NewId, pGeometry) {}

template <class TElementData>
TwoFluidNavierStokes<TElementData>::TwoFluidNavierStokes(
    IndexType NewId, GeometryType::Pointer pGeometry, Properties::Pointer pProperties)
    : FluidElement<TElementData>(NewId, pGeometry, pProperties) {}

template <class TElementData>
TwoFluidNavierStokes<TElementData>::~TwoFluidNavierStokes() {}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template <class TElementData>
Element::Pointer TwoFluidNavierStokes<TElementData>::Create(
    IndexType NewId,
    NodesArrayType const &ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<TwoFluidNavierStokes>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}

template <class TElementData>
Element::Pointer TwoFluidNavierStokes<TElementData>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<TwoFluidNavierStokes>(NewId, pGeom, pProperties);
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::CalculateLocalSystem(
    MatrixType &rLeftHandSideMatrix,
    VectorType &rRightHandSideVector,
    const ProcessInfo &rCurrentProcessInfo)
{
    // Resize and intialize output
    if (rLeftHandSideMatrix.size1() != LocalSize)
        rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);

    if (rRightHandSideVector.size() != LocalSize)
        rRightHandSideVector.resize(LocalSize, false);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);
    noalias(rRightHandSideVector) = ZeroVector(LocalSize);

    if (TElementData::ElementManagesTimeIntegration){
        TElementData data;
        data.Initialize(*this, rCurrentProcessInfo);

        if (data.IsCut()){
            GeometryType::Pointer p_geom = this->pGetGeometry();
            Matrix shape_functions_pos, shape_functions_neg;
            Matrix shape_functions_enr_pos, shape_functions_enr_neg;
            GeometryType::ShapeFunctionsGradientsType shape_derivatives_pos, shape_derivatives_neg;
            GeometryType::ShapeFunctionsGradientsType shape_derivatives_enr_pos, shape_derivatives_enr_neg;

            ModifiedShapeFunctions::Pointer p_modified_sh_func = pGetModifiedShapeFunctionsUtility(p_geom, data.Distance);

            ComputeSplitting(
                data,
                shape_functions_pos,
                shape_functions_neg,
                shape_functions_enr_pos,
                shape_functions_enr_neg,
                shape_derivatives_pos,
                shape_derivatives_neg,
                shape_derivatives_enr_pos,
                shape_derivatives_enr_neg,
                p_modified_sh_func);

            if (data.NumberOfDivisions == 1){
                // Cases exist when the element is not subdivided due to the characteristics of the provided distance
                // In this cases the element is treated as AIR or FLUID depending on the side
                Vector gauss_weights;
                Matrix shape_functions;
                ShapeFunctionDerivativesArrayType shape_derivatives;
                this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
                const unsigned int number_of_gauss_points = gauss_weights.size();
                array_1d<double, NumNodes> Ncenter;
                for (unsigned int i = 0; i < NumNodes; ++i){
                    Ncenter[i] = 1.0/NumNodes;
                }
                for (unsigned int g = 0; g < number_of_gauss_points; ++g){
                    UpdateIntegrationPointData(
                        data,
                        g,
                        gauss_weights[g],
                        row(shape_functions, g),
                        shape_derivatives[g]);
                    this->AddTimeIntegratedSystem(data, rLeftHandSideMatrix, rRightHandSideVector);
                }
            } else {
                MatrixType Vtot = ZeroMatrix(NumNodes * (Dim + 1), NumNodes);
                MatrixType Htot = ZeroMatrix(NumNodes, NumNodes * (Dim + 1));
                MatrixType Kee_tot = ZeroMatrix(NumNodes, NumNodes);
                VectorType rhs_ee_tot = ZeroVector(NumNodes);

                for (unsigned int g_pos = 0; g_pos < data.w_gauss_pos_side.size(); ++g_pos){
                    UpdateIntegrationPointData(
                        data,
                        g_pos,
                        data.w_gauss_pos_side[g_pos],
                        row(shape_functions_pos, g_pos),
                        shape_derivatives_pos[g_pos],
                        row(shape_functions_enr_pos, g_pos),
                        shape_derivatives_enr_pos[g_pos]);

                    this->AddTimeIntegratedSystem(data, rLeftHandSideMatrix, rRightHandSideVector);
                    this->ComputeGaussPointEnrichmentContributions(data, Vtot, Htot, Kee_tot, rhs_ee_tot);
                }

                for (unsigned int g_neg = 0; g_neg < data.w_gauss_neg_side.size(); ++g_neg){
                    UpdateIntegrationPointData(
                        data,
                        g_neg,
                        data.w_gauss_neg_side[g_neg],
                        row(shape_functions_neg, g_neg),
                        shape_derivatives_neg[g_neg],
                        row(shape_functions_enr_neg, g_neg),
                        shape_derivatives_enr_neg[g_neg]);
                    this->AddTimeIntegratedSystem(data, rLeftHandSideMatrix, rRightHandSideVector);
                    this->ComputeGaussPointEnrichmentContributions(data, Vtot, Htot, Kee_tot, rhs_ee_tot);
                }

                Matrix int_shape_function, int_shape_function_enr_neg, int_shape_function_enr_pos;
                GeometryType::ShapeFunctionsGradientsType int_shape_derivatives;
                Vector int_gauss_pts_weights;
                std::vector< array_1d<double,3> > int_normals_neg;

                // Base momentum correction and momentum/mass correction are incompatible
                KRATOS_ERROR_IF(rCurrentProcessInfo[MOMENTUM_CORRECTION] && rCurrentProcessInfo[ENERGY_PRESERVING_MOMENTUM_TERM]) << "Base momentum correction and momentum/mass correction are incompatible." << std::endl;

                if (rCurrentProcessInfo[SURFACE_TENSION] || rCurrentProcessInfo[MOMENTUM_CORRECTION] || rCurrentProcessInfo[ENERGY_PRESERVING_MOMENTUM_TERM])
                {
                    ComputeSplitInterface(
                        data,
                        int_shape_function,
                        int_shape_function_enr_pos,
                        int_shape_function_enr_neg,
                        int_shape_derivatives,
                        int_gauss_pts_weights,
                        int_normals_neg,
                        p_modified_sh_func);
                }
                if (rCurrentProcessInfo[ENERGY_PRESERVING_MOMENTUM_TERM]){

                    BoundedMatrix<double, LocalSize, LocalSize> lhs_acc_correction = ZeroMatrix(LocalSize, LocalSize);
                    double positive_density = 0.0;
                    double negative_density = 0.0;
                    const auto &r_geom = this->GetGeometry();
                    for (unsigned int intgp = 0; intgp < int_gauss_pts_weights.size(); ++intgp)
                    {
                        double u_distance = 0.0;
                        array_1d<double, 3> u_field = ZeroVector(3);
                        for (unsigned int i = 0; i < NumNodes; ++i)
                        {
                            noalias(u_field) += int_shape_function(intgp, i) * r_geom[i].FastGetSolutionStepValue(VELOCITY, 1);
                            u_distance += int_shape_function(intgp, i) * r_geom[i].FastGetSolutionStepValue(DISTANCE, 1);

                            if (data.Distance[i] > 0.0)
                            {
                                positive_density = data.NodalDensity[i];
                            }
                            else
                            {
                                negative_density = data.NodalDensity[i];
                            }
                        }

                        u_distance /= data.DeltaTime;
                        const array_1d<double, 3> &r_n = int_normals_neg[intgp];
                        double u_prime = u_distance;
                        u_prime -= inner_prod(u_field, r_n);

                        for (unsigned int i = 0; i < BlockSize; ++i)
                        {
                            for (unsigned int j = 0; j < BlockSize; ++j)
                            {
                                for (unsigned int dim = 0; dim < Dim; ++dim)
                                {
                                    lhs_acc_correction(i * (BlockSize) + dim, j * (BlockSize) + dim) +=
                                    int_shape_function(intgp, i) * int_shape_function(intgp, j) * u_prime * int_gauss_pts_weights(intgp);
                                }
                            }
                        }
                    }

                    lhs_acc_correction *= (positive_density-negative_density);
                    noalias(rLeftHandSideMatrix) += lhs_acc_correction;

                    Kratos::array_1d<double, LocalSize> tempU; // Unknowns vector containing only velocity components
                    for (unsigned int i = 0; i < NumNodes; ++i)
                    {
                        for (unsigned int dimi = 0; dimi < Dim; ++dimi)
                        {
                            tempU[i * (BlockSize) + dimi] = data.Velocity(i, dimi);
                        }
                    }
                    noalias(rRightHandSideVector) -= prod(lhs_acc_correction, tempU);
                }

                if (rCurrentProcessInfo[MOMENTUM_CORRECTION]){
                    BoundedMatrix<double, LocalSize, LocalSize> lhs_acc_correction = ZeroMatrix(LocalSize,LocalSize);

                    double positive_density = 0.0;
                    double negative_density = 0.0;

                    const auto& r_geom = this->GetGeometry();

                    for (unsigned int intgp = 0; intgp < int_gauss_pts_weights.size(); ++intgp){
                        double u_dot_n = 0.0;
                        for (unsigned int i = 0; i < NumNodes; ++i){
                            u_dot_n += int_shape_function(intgp,i)*r_geom[i].GetValue(DISTANCE_CORRECTION);

                            if (data.Distance[i] > 0.0){
                                positive_density = data.NodalDensity[i];
                            } else {
                                negative_density = data.NodalDensity[i];
                            }
                        }

                        u_dot_n /= data.DeltaTime;

                        for (unsigned int i = 0; i < NumNodes; ++i){
                            for (unsigned int j = 0; j < NumNodes; ++j){
                                for (unsigned int dim = 0; dim < NumNodes-1; ++dim){
                                    lhs_acc_correction( i*(NumNodes) + dim, j*(NumNodes) + dim) +=
                                        int_shape_function(intgp,i)*int_shape_function(intgp,j)*u_dot_n*int_gauss_pts_weights(intgp);
                                }
                            }
                        }
                    }

                    lhs_acc_correction = (negative_density - positive_density)*lhs_acc_correction;
                    noalias(rLeftHandSideMatrix) += lhs_acc_correction;

                    Kratos::array_1d<double, LocalSize> tempU; // Unknowns vector containing only velocity components
                    for (unsigned int i = 0; i < NumNodes; ++i){
                        for (unsigned int dimi = 0; dimi < Dim; ++dimi){
                            tempU[i*(Dim+1) + dimi] = data.Velocity(i,dimi);
                        }
                    }
                    noalias(rRightHandSideVector) -= prod(lhs_acc_correction,tempU);
                }

                if (rCurrentProcessInfo[SURFACE_TENSION]){

                    AddSurfaceTensionContribution(
                        data,
                        int_shape_function,
                        int_shape_function_enr_pos,
                        int_shape_function_enr_neg,
                        int_shape_derivatives,
                        int_gauss_pts_weights,
                        int_normals_neg,
                        rLeftHandSideMatrix,
                        rRightHandSideVector,
                        Htot,
                        Vtot,
                        Kee_tot,
                        rhs_ee_tot
                    );

                } else{
                    // Without pressure gradient stabilization, volume ratio is checked during condensation
                    // Also, without surface tension, zero pressure difference is penalized
                    CondenseEnrichmentWithContinuity(data, rLeftHandSideMatrix, rRightHandSideVector, Htot, Vtot, Kee_tot, rhs_ee_tot);
                }

            }
        } else {
            //Get Shape function data
            Vector gauss_weights;
            Matrix shape_functions;
            ShapeFunctionDerivativesArrayType shape_derivatives;
            this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
            const unsigned int number_of_gauss_points = gauss_weights.size();
            // Iterate over integration points to evaluate local contribution
            for (unsigned int g = 0; g < number_of_gauss_points; ++g){
                UpdateIntegrationPointData(data, g, gauss_weights[g], row(shape_functions, g), shape_derivatives[g]);
                this->AddTimeIntegratedSystem(data, rLeftHandSideMatrix, rRightHandSideVector);
            }
        }
    } else{
        KRATOS_ERROR << "TwoFluidNavierStokes is supposed to manage time integration." << std::endl;
    }
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::CalculateRightHandSide(
    VectorType &rRightHandSideVector,
    const ProcessInfo &rCurrentProcessInfo)
{
    MatrixType tmp;
    CalculateLocalSystem(tmp, rRightHandSideVector, rCurrentProcessInfo);
}
///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Inquiry

template <class TElementData>
int TwoFluidNavierStokes<TElementData>::Check(const ProcessInfo &rCurrentProcessInfo) const
{
    KRATOS_TRY;
    int out = FluidElement<TElementData>::Check(rCurrentProcessInfo);
    KRATOS_ERROR_IF_NOT(out == 0)
        << "Error in base class Check for Element " << this->Info() << std::endl
        << "Error code is " << out << std::endl;

    return 0;

    KRATOS_CATCH("");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public I/O

template <class TElementData>
const Parameters TwoFluidNavierStokes<TElementData>::GetSpecifications() const
{
    const Parameters specifications = Parameters(R"({
        "time_integration"           : ["implicit"],
        "framework"                  : "ale",
        "symmetric_lhs"              : false,
        "positive_definite_lhs"      : true,
        "output"                     : {
            "gauss_point"            : [],
            "nodal_historical"       : ["VELOCITY","PRESSURE"],
            "nodal_non_historical"   : [],
            "entity"                 : []
        },
        "required_variables"         : ["DISTANCE","VELOCITY","PRESSURE","MESH_VELOCITY","DENSITY","DYNAMIC_VISCOSITY"],
        "required_dofs"              : [],
        "flags_used"                 : [],
        "compatible_geometries"      : ["Triangle2D3","Tetrahedra3D4"],
        "element_integrates_in_time" : true,
        "compatible_constitutive_laws": {
            "type"        : ["NewtonianTwoFluid2DLaw","NewtonianTwoFluid3DLaw"],
            "dimension"   : ["2D","3D"],
            "strain_size" : [3,6]
        },
        "required_polynomial_degree_of_geometry" : 1,
        "documentation"   :
            "This element implements Navier-Stokes biphasic fluid-air formulation with a levelset-based interface representation with Variational MultiScales (VMS) stabilization. Note that any viscous behavior can be used for the fluid phase through a constitutive law. The air phase is assumed to be Newtonian. Surface tension contribution can be accounted for by setting the SURFACE_TENSION variable to true in the ProcessInfo container.
    })");

    if (Dim == 2) {
        std::vector<std::string> dofs_2d({"VELOCITY_X","VELOCITY_Y","PRESSURE"});
        specifications["required_dofs"].SetStringArray(dofs_2d);
    } else {
        std::vector<std::string> dofs_3d({"VELOCITY_X","VELOCITY_Y","VELOCITY_Z","PRESSURE"});
        specifications["required_dofs"].SetStringArray(dofs_3d);
    }

    return specifications;
}

template <class TElementData>
std::string TwoFluidNavierStokes<TElementData>::Info() const
{
    std::stringstream buffer;
    buffer << "TwoFluidNavierStokes" << Dim << "D" << NumNodes << "N #" << this->Id();
    return buffer.str();
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::PrintInfo(
    std::ostream &rOStream) const
{
    rOStream << this->Info() << std::endl;

    if (this->GetConstitutiveLaw() != nullptr){
        rOStream << "with constitutive law " << std::endl;
        this->GetConstitutiveLaw()->PrintInfo(rOStream);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected operations

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::AddTimeIntegratedSystem(
    TElementData &rData,
    MatrixType &rLHS,
    VectorType &rRHS)
{
    this->ComputeGaussPointLHSContribution(rData, rLHS);
    this->ComputeGaussPointRHSContribution(rData, rRHS);
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::AddTimeIntegratedLHS(
    TElementData &rData,
    MatrixType &rLHS)
{
    this->ComputeGaussPointLHSContribution(rData, rLHS);
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::AddTimeIntegratedRHS(
    TElementData &rData,
    VectorType &rRHS)
{
    this->ComputeGaussPointRHSContribution(rData, rRHS);
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::UpdateIntegrationPointData(
    TElementData& rData,
    unsigned int IntegrationPointIndex,
    double Weight,
    const typename TElementData::MatrixRowType& rN,
    const typename TElementData::ShapeDerivativesType& rDN_DX) const
{
    rData.UpdateGeometryValues(IntegrationPointIndex, Weight, rN, rDN_DX);
    const double d_gauss = inner_prod(rData.Distance, rN);
    if (d_gauss > 0.0)
        rData.CalculateAirMaterialResponse();
    else
        this->CalculateMaterialResponse(rData);
    rData.ComputeDarcyTerm();
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::UpdateIntegrationPointData(
    TElementData& rData,
    unsigned int IntegrationPointIndex,
    double Weight,
    const typename TElementData::MatrixRowType& rN,
    const typename TElementData::ShapeDerivativesType& rDN_DX,
    const typename TElementData::MatrixRowType& rNenr,
    const typename TElementData::ShapeDerivativesType& rDN_DXenr) const
{
    rData.UpdateGeometryValues(IntegrationPointIndex,Weight,rN,rDN_DX,rNenr,rDN_DXenr);
    const double d_gauss = inner_prod(rData.Distance, rN);
    if (d_gauss > 0.0)
        rData.CalculateAirMaterialResponse();
    else
        this->CalculateMaterialResponse(rData);
    rData.ComputeDarcyTerm();
}

template <>
void TwoFluidNavierStokes<TwoFluidNavierStokesData<2, 3>>::ComputeGaussPointLHSContribution(
    TwoFluidNavierStokesData<2, 3> &rData,
    MatrixType &rLHS)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;

    const double dyn_tau = rData.DynamicTau;
    const double K_darcy = rData.DarcyTerm;

    const auto vconv = rData.Velocity - rData.MeshVelocity;
    const double volume_error_ratio = rData.VolumeError;
    // Get constitutive matrix
    const Matrix &C = rData.C;

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    auto &lhs = rData.lhs;

    const double clhs0 = C(0,0)*DN(0,0) + C(0,2)*DN(0,1);
const double clhs1 = C(0,2)*DN(0,0);
const double clhs2 = C(2,2)*DN(0,1) + clhs1;
const double clhs3 = pow(DN(0,0), 2);
const double clhs4 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double clhs5 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double clhs6 = rho*stab_c2*sqrt(pow(clhs4, 2) + pow(clhs5, 2));
const double clhs7 = clhs6*h/stab_c1 + mu;
const double clhs8 = pow(N[0], 2);
const double clhs9 = rho*(DN(0,0)*clhs4 + DN(0,1)*clhs5);
const double clhs10 = clhs8*rho;
const double clhs11 = K_darcy*N[0];
const double clhs12 = N[0]*rho;
const double clhs13 = bdf0*clhs12;
const double clhs14 = clhs12*volume_error_ratio;
const double clhs15 = clhs11 + clhs13 + clhs14 + clhs9;
const double clhs16 = 1.0/(K_darcy + clhs6/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs17 = 1.0*clhs9;
const double clhs18 = clhs16*clhs17;
const double clhs19 = 1.0*clhs11;
const double clhs20 = clhs16*clhs19;
const double clhs21 = 1.0*clhs16;
const double clhs22 = clhs15*clhs21;
const double clhs23 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double clhs24 = clhs12*clhs23;
const double clhs25 = K_darcy*clhs8 + N[0]*clhs9 + bdf0*clhs10 + clhs10*volume_error_ratio + clhs15*clhs18 - clhs15*clhs20 + clhs22*clhs24;
const double clhs26 = C(0,1)*DN(0,1) + clhs1;
const double clhs27 = C(1,2)*DN(0,1);
const double clhs28 = C(2,2)*DN(0,0) + clhs27;
const double clhs29 = DN(0,0)*clhs7;
const double clhs30 = DN(0,1)*clhs29;
const double clhs31 = clhs21*clhs23;
const double clhs32 = -N[0] + clhs12*clhs31 + clhs18 - clhs20;
const double clhs33 = C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double clhs34 = C(0,2)*DN(1,0);
const double clhs35 = C(2,2)*DN(1,1) + clhs34;
const double clhs36 = DN(0,0)*DN(1,0);
const double clhs37 = N[1]*clhs11 + N[1]*clhs13 + N[1]*clhs14;
const double clhs38 = clhs36*clhs7 + clhs37;
const double clhs39 = rho*(DN(1,0)*clhs4 + DN(1,1)*clhs5);
const double clhs40 = K_darcy*N[1];
const double clhs41 = N[1]*rho;
const double clhs42 = bdf0*clhs41;
const double clhs43 = clhs41*volume_error_ratio;
const double clhs44 = clhs39 + clhs40 + clhs42 + clhs43;
const double clhs45 = clhs21*clhs44;
const double clhs46 = N[0]*clhs39 + clhs18*clhs44 - clhs20*clhs44 + clhs24*clhs45;
const double clhs47 = C(0,1)*DN(1,1) + clhs34;
const double clhs48 = C(1,2)*DN(1,1);
const double clhs49 = C(2,2)*DN(1,0) + clhs48;
const double clhs50 = DN(1,1)*clhs29;
const double clhs51 = DN(0,0)*N[1];
const double clhs52 = DN(1,0)*N[0];
const double clhs53 = clhs31*rho;
const double clhs54 = C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double clhs55 = C(0,2)*DN(2,0);
const double clhs56 = C(2,2)*DN(2,1) + clhs55;
const double clhs57 = DN(0,0)*DN(2,0);
const double clhs58 = N[2]*clhs11 + N[2]*clhs13 + N[2]*clhs14;
const double clhs59 = clhs57*clhs7 + clhs58;
const double clhs60 = rho*(DN(2,0)*clhs4 + DN(2,1)*clhs5);
const double clhs61 = K_darcy*N[2];
const double clhs62 = N[2]*rho;
const double clhs63 = bdf0*clhs62;
const double clhs64 = clhs62*volume_error_ratio;
const double clhs65 = clhs60 + clhs61 + clhs63 + clhs64;
const double clhs66 = clhs21*clhs65;
const double clhs67 = N[0]*clhs60 + clhs18*clhs65 - clhs20*clhs65 + clhs24*clhs66;
const double clhs68 = C(0,1)*DN(2,1) + clhs55;
const double clhs69 = C(1,2)*DN(2,1);
const double clhs70 = C(2,2)*DN(2,0) + clhs69;
const double clhs71 = DN(2,1)*clhs29;
const double clhs72 = DN(0,0)*N[2];
const double clhs73 = DN(2,0)*N[0];
const double clhs74 = C(0,1)*DN(0,0) + clhs27;
const double clhs75 = C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double clhs76 = pow(DN(0,1), 2);
const double clhs77 = C(0,1)*DN(1,0) + clhs48;
const double clhs78 = DN(0,1)*clhs7;
const double clhs79 = DN(1,0)*clhs78;
const double clhs80 = C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double clhs81 = DN(0,1)*DN(1,1);
const double clhs82 = clhs37 + clhs7*clhs81;
const double clhs83 = DN(0,1)*N[1];
const double clhs84 = DN(1,1)*N[0];
const double clhs85 = C(0,1)*DN(2,0) + clhs69;
const double clhs86 = DN(2,0)*clhs78;
const double clhs87 = C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double clhs88 = DN(0,1)*DN(2,1);
const double clhs89 = clhs58 + clhs7*clhs88;
const double clhs90 = DN(0,1)*N[2];
const double clhs91 = DN(2,1)*N[0];
const double clhs92 = N[0] + clhs16*(1.0*clhs13 + 1.0*clhs14 + clhs17 + clhs19);
const double clhs93 = clhs21*(clhs36 + clhs81);
const double clhs94 = clhs21*(clhs57 + clhs88);
const double clhs95 = clhs21*clhs39;
const double clhs96 = clhs21*clhs40;
const double clhs97 = clhs23*clhs41;
const double clhs98 = N[1]*clhs9 + clhs15*clhs95 - clhs15*clhs96 + clhs22*clhs97;
const double clhs99 = pow(DN(1,0), 2);
const double clhs100 = pow(N[1], 2);
const double clhs101 = clhs100*rho;
const double clhs102 = K_darcy*clhs100 + N[1]*clhs39 + bdf0*clhs101 + clhs101*volume_error_ratio + clhs39*clhs45 - clhs40*clhs45 + clhs45*clhs97;
const double clhs103 = DN(1,0)*clhs7;
const double clhs104 = DN(1,1)*clhs103;
const double clhs105 = -N[1] + clhs31*clhs41 + clhs95 - clhs96;
const double clhs106 = DN(1,0)*DN(2,0);
const double clhs107 = N[2]*clhs40 + N[2]*clhs42 + N[2]*clhs43;
const double clhs108 = clhs106*clhs7 + clhs107;
const double clhs109 = N[1]*clhs60 + clhs39*clhs66 - clhs40*clhs66 + clhs66*clhs97;
const double clhs110 = DN(2,1)*clhs103;
const double clhs111 = DN(1,0)*N[2];
const double clhs112 = DN(2,0)*N[1];
const double clhs113 = pow(DN(1,1), 2);
const double clhs114 = DN(2,0)*clhs7;
const double clhs115 = DN(1,1)*clhs114;
const double clhs116 = DN(1,1)*DN(2,1);
const double clhs117 = clhs107 + clhs116*clhs7;
const double clhs118 = DN(1,1)*N[2];
const double clhs119 = DN(2,1)*N[1];
const double clhs120 = N[1] + clhs16*(1.0*clhs39 + 1.0*clhs40 + 1.0*clhs42 + 1.0*clhs43);
const double clhs121 = clhs21*(clhs106 + clhs116);
const double clhs122 = clhs23*clhs62;
const double clhs123 = N[2]*clhs9 + clhs122*clhs22 + clhs22*clhs60 - clhs22*clhs61;
const double clhs124 = clhs21*clhs61;
const double clhs125 = clhs21*clhs60;
const double clhs126 = N[2]*clhs39 + clhs122*clhs45 + clhs45*clhs60 - clhs45*clhs61;
const double clhs127 = pow(DN(2,0), 2);
const double clhs128 = pow(N[2], 2);
const double clhs129 = clhs128*rho;
const double clhs130 = K_darcy*clhs128 + N[2]*clhs60 + bdf0*clhs129 + clhs122*clhs66 + clhs129*volume_error_ratio + clhs60*clhs66 - clhs61*clhs66;
const double clhs131 = DN(2,1)*clhs114;
const double clhs132 = -N[2] - clhs124 + clhs125 + clhs31*clhs62;
const double clhs133 = pow(DN(2,1), 2);
const double clhs134 = N[2] + clhs16*(1.0*clhs60 + 1.0*clhs61 + 1.0*clhs63 + 1.0*clhs64);
lhs(0,0)=DN(0,0)*clhs0 + DN(0,1)*clhs2 + clhs25 + clhs3*clhs7;
lhs(0,1)=DN(0,0)*clhs26 + DN(0,1)*clhs28 + clhs30;
lhs(0,2)=DN(0,0)*clhs32;
lhs(0,3)=DN(0,0)*clhs33 + DN(0,1)*clhs35 + clhs38 + clhs46;
lhs(0,4)=DN(0,0)*clhs47 + DN(0,1)*clhs49 + clhs50;
lhs(0,5)=DN(1,0)*clhs18 - DN(1,0)*clhs20 - clhs51 + clhs52*clhs53;
lhs(0,6)=DN(0,0)*clhs54 + DN(0,1)*clhs56 + clhs59 + clhs67;
lhs(0,7)=DN(0,0)*clhs68 + DN(0,1)*clhs70 + clhs71;
lhs(0,8)=DN(2,0)*clhs18 - DN(2,0)*clhs20 + clhs53*clhs73 - clhs72;
lhs(1,0)=DN(0,0)*clhs2 + DN(0,1)*clhs74 + clhs30;
lhs(1,1)=DN(0,0)*clhs28 + DN(0,1)*clhs75 + clhs25 + clhs7*clhs76;
lhs(1,2)=DN(0,1)*clhs32;
lhs(1,3)=DN(0,0)*clhs35 + DN(0,1)*clhs77 + clhs79;
lhs(1,4)=DN(0,0)*clhs49 + DN(0,1)*clhs80 + clhs46 + clhs82;
lhs(1,5)=DN(1,1)*clhs18 - DN(1,1)*clhs20 + clhs53*clhs84 - clhs83;
lhs(1,6)=DN(0,0)*clhs56 + DN(0,1)*clhs85 + clhs86;
lhs(1,7)=DN(0,0)*clhs70 + DN(0,1)*clhs87 + clhs67 + clhs89;
lhs(1,8)=DN(2,1)*clhs18 - DN(2,1)*clhs20 + clhs53*clhs91 - clhs90;
lhs(2,0)=DN(0,0)*clhs92;
lhs(2,1)=DN(0,1)*clhs92;
lhs(2,2)=clhs21*(clhs3 + clhs76);
lhs(2,3)=DN(0,0)*clhs45 + clhs52;
lhs(2,4)=DN(0,1)*clhs45 + clhs84;
lhs(2,5)=clhs93;
lhs(2,6)=DN(0,0)*clhs66 + clhs73;
lhs(2,7)=DN(0,1)*clhs66 + clhs91;
lhs(2,8)=clhs94;
lhs(3,0)=DN(1,0)*clhs0 + DN(1,1)*clhs2 + clhs38 + clhs98;
lhs(3,1)=DN(1,0)*clhs26 + DN(1,1)*clhs28 + clhs79;
lhs(3,2)=DN(0,0)*clhs95 - DN(0,0)*clhs96 + clhs51*clhs53 - clhs52;
lhs(3,3)=DN(1,0)*clhs33 + DN(1,1)*clhs35 + clhs102 + clhs7*clhs99;
lhs(3,4)=DN(1,0)*clhs47 + DN(1,1)*clhs49 + clhs104;
lhs(3,5)=DN(1,0)*clhs105;
lhs(3,6)=DN(1,0)*clhs54 + DN(1,1)*clhs56 + clhs108 + clhs109;
lhs(3,7)=DN(1,0)*clhs68 + DN(1,1)*clhs70 + clhs110;
lhs(3,8)=DN(2,0)*clhs95 - DN(2,0)*clhs96 - clhs111 + clhs112*clhs53;
lhs(4,0)=DN(1,0)*clhs2 + DN(1,1)*clhs74 + clhs50;
lhs(4,1)=DN(1,0)*clhs28 + DN(1,1)*clhs75 + clhs82 + clhs98;
lhs(4,2)=DN(0,1)*clhs95 - DN(0,1)*clhs96 + clhs53*clhs83 - clhs84;
lhs(4,3)=DN(1,0)*clhs35 + DN(1,1)*clhs77 + clhs104;
lhs(4,4)=DN(1,0)*clhs49 + DN(1,1)*clhs80 + clhs102 + clhs113*clhs7;
lhs(4,5)=DN(1,1)*clhs105;
lhs(4,6)=DN(1,0)*clhs56 + DN(1,1)*clhs85 + clhs115;
lhs(4,7)=DN(1,0)*clhs70 + DN(1,1)*clhs87 + clhs109 + clhs117;
lhs(4,8)=DN(2,1)*clhs95 - DN(2,1)*clhs96 - clhs118 + clhs119*clhs53;
lhs(5,0)=DN(1,0)*clhs22 + clhs51;
lhs(5,1)=DN(1,1)*clhs22 + clhs83;
lhs(5,2)=clhs93;
lhs(5,3)=DN(1,0)*clhs120;
lhs(5,4)=DN(1,1)*clhs120;
lhs(5,5)=clhs21*(clhs113 + clhs99);
lhs(5,6)=DN(1,0)*clhs66 + clhs112;
lhs(5,7)=DN(1,1)*clhs66 + clhs119;
lhs(5,8)=clhs121;
lhs(6,0)=DN(2,0)*clhs0 + DN(2,1)*clhs2 + clhs123 + clhs59;
lhs(6,1)=DN(2,0)*clhs26 + DN(2,1)*clhs28 + clhs86;
lhs(6,2)=-DN(0,0)*clhs124 + DN(0,0)*clhs125 + clhs53*clhs72 - clhs73;
lhs(6,3)=DN(2,0)*clhs33 + DN(2,1)*clhs35 + clhs108 + clhs126;
lhs(6,4)=DN(2,0)*clhs47 + DN(2,1)*clhs49 + clhs115;
lhs(6,5)=-DN(1,0)*clhs124 + DN(1,0)*clhs125 + clhs111*clhs53 - clhs112;
lhs(6,6)=DN(2,0)*clhs54 + DN(2,1)*clhs56 + clhs127*clhs7 + clhs130;
lhs(6,7)=DN(2,0)*clhs68 + DN(2,1)*clhs70 + clhs131;
lhs(6,8)=DN(2,0)*clhs132;
lhs(7,0)=DN(2,0)*clhs2 + DN(2,1)*clhs74 + clhs71;
lhs(7,1)=DN(2,0)*clhs28 + DN(2,1)*clhs75 + clhs123 + clhs89;
lhs(7,2)=-DN(0,1)*clhs124 + DN(0,1)*clhs125 + clhs53*clhs90 - clhs91;
lhs(7,3)=DN(2,0)*clhs35 + DN(2,1)*clhs77 + clhs110;
lhs(7,4)=DN(2,0)*clhs49 + DN(2,1)*clhs80 + clhs117 + clhs126;
lhs(7,5)=-DN(1,1)*clhs124 + DN(1,1)*clhs125 + clhs118*clhs53 - clhs119;
lhs(7,6)=DN(2,0)*clhs56 + DN(2,1)*clhs85 + clhs131;
lhs(7,7)=DN(2,0)*clhs70 + DN(2,1)*clhs87 + clhs130 + clhs133*clhs7;
lhs(7,8)=DN(2,1)*clhs132;
lhs(8,0)=DN(2,0)*clhs22 + clhs72;
lhs(8,1)=DN(2,1)*clhs22 + clhs90;
lhs(8,2)=clhs94;
lhs(8,3)=DN(2,0)*clhs45 + clhs111;
lhs(8,4)=DN(2,1)*clhs45 + clhs118;
lhs(8,5)=clhs121;
lhs(8,6)=DN(2,0)*clhs134;
lhs(8,7)=DN(2,1)*clhs134;
lhs(8,8)=clhs21*(clhs127 + clhs133);


    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;
}

template <>
void TwoFluidNavierStokes<TwoFluidNavierStokesData<3, 4>>::ComputeGaussPointLHSContribution(
    TwoFluidNavierStokesData<3, 4> &rData,
    MatrixType &rLHS)
{

    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;
    const double K_darcy = rData.DarcyTerm;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;

    const double dyn_tau = rData.DynamicTau;

    const auto vconv = rData.Velocity - rData.MeshVelocity;
    const double volume_error_ratio = rData.VolumeError;
    // Get constitutive matrix
    const Matrix &C = rData.C;

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    auto &lhs = rData.lhs;

    const double clhs0 = C(0,0)*DN(0,0) + C(0,3)*DN(0,1) + C(0,5)*DN(0,2);
const double clhs1 = C(0,3)*DN(0,0);
const double clhs2 = C(3,3)*DN(0,1) + C(3,5)*DN(0,2) + clhs1;
const double clhs3 = C(0,5)*DN(0,0);
const double clhs4 = C(3,5)*DN(0,1) + C(5,5)*DN(0,2) + clhs3;
const double clhs5 = pow(DN(0,0), 2);
const double clhs6 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double clhs7 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double clhs8 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double clhs9 = rho*stab_c2*sqrt(pow(clhs6, 2) + pow(clhs7, 2) + pow(clhs8, 2));
const double clhs10 = clhs9*h/stab_c1 + mu;
const double clhs11 = pow(N[0], 2);
const double clhs12 = rho*(DN(0,0)*clhs6 + DN(0,1)*clhs7 + DN(0,2)*clhs8);
const double clhs13 = clhs11*rho;
const double clhs14 = K_darcy*N[0];
const double clhs15 = N[0]*rho;
const double clhs16 = bdf0*clhs15;
const double clhs17 = clhs15*volume_error_ratio;
const double clhs18 = clhs12 + clhs14 + clhs16 + clhs17;
const double clhs19 = 1.0/(K_darcy + clhs9/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs20 = 1.0*clhs12;
const double clhs21 = clhs19*clhs20;
const double clhs22 = 1.0*clhs14;
const double clhs23 = clhs19*clhs22;
const double clhs24 = 1.0*clhs19;
const double clhs25 = clhs18*clhs24;
const double clhs26 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double clhs27 = clhs15*clhs26;
const double clhs28 = K_darcy*clhs11 + N[0]*clhs12 + bdf0*clhs13 + clhs13*volume_error_ratio + clhs18*clhs21 - clhs18*clhs23 + clhs25*clhs27;
const double clhs29 = C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + clhs1;
const double clhs30 = C(1,3)*DN(0,1);
const double clhs31 = C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + clhs30;
const double clhs32 = C(3,5)*DN(0,0);
const double clhs33 = C(4,5)*DN(0,2);
const double clhs34 = C(1,5)*DN(0,1) + clhs32 + clhs33;
const double clhs35 = DN(0,0)*clhs10;
const double clhs36 = DN(0,1)*clhs35;
const double clhs37 = C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + clhs3;
const double clhs38 = C(3,4)*DN(0,1);
const double clhs39 = C(2,3)*DN(0,2) + clhs32 + clhs38;
const double clhs40 = C(2,5)*DN(0,2);
const double clhs41 = C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + clhs40;
const double clhs42 = DN(0,2)*clhs35;
const double clhs43 = clhs24*clhs26;
const double clhs44 = -N[0] + clhs15*clhs43 + clhs21 - clhs23;
const double clhs45 = C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double clhs46 = C(0,3)*DN(1,0);
const double clhs47 = C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + clhs46;
const double clhs48 = C(0,5)*DN(1,0);
const double clhs49 = C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + clhs48;
const double clhs50 = DN(0,0)*DN(1,0);
const double clhs51 = N[1]*clhs14 + N[1]*clhs16 + N[1]*clhs17;
const double clhs52 = clhs10*clhs50 + clhs51;
const double clhs53 = rho*(DN(1,0)*clhs6 + DN(1,1)*clhs7 + DN(1,2)*clhs8);
const double clhs54 = K_darcy*N[1];
const double clhs55 = N[1]*rho;
const double clhs56 = bdf0*clhs55;
const double clhs57 = clhs55*volume_error_ratio;
const double clhs58 = clhs53 + clhs54 + clhs56 + clhs57;
const double clhs59 = clhs24*clhs58;
const double clhs60 = N[0]*clhs53 + clhs21*clhs58 - clhs23*clhs58 + clhs27*clhs59;
const double clhs61 = C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + clhs46;
const double clhs62 = C(1,3)*DN(1,1);
const double clhs63 = C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + clhs62;
const double clhs64 = C(3,5)*DN(1,0);
const double clhs65 = C(4,5)*DN(1,2);
const double clhs66 = C(1,5)*DN(1,1) + clhs64 + clhs65;
const double clhs67 = DN(1,1)*clhs35;
const double clhs68 = C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + clhs48;
const double clhs69 = C(3,4)*DN(1,1);
const double clhs70 = C(2,3)*DN(1,2) + clhs64 + clhs69;
const double clhs71 = C(2,5)*DN(1,2);
const double clhs72 = C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + clhs71;
const double clhs73 = DN(1,2)*clhs35;
const double clhs74 = DN(0,0)*N[1];
const double clhs75 = DN(1,0)*N[0];
const double clhs76 = clhs43*rho;
const double clhs77 = C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double clhs78 = C(0,3)*DN(2,0);
const double clhs79 = C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + clhs78;
const double clhs80 = C(0,5)*DN(2,0);
const double clhs81 = C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + clhs80;
const double clhs82 = DN(0,0)*DN(2,0);
const double clhs83 = N[2]*clhs14 + N[2]*clhs16 + N[2]*clhs17;
const double clhs84 = clhs10*clhs82 + clhs83;
const double clhs85 = rho*(DN(2,0)*clhs6 + DN(2,1)*clhs7 + DN(2,2)*clhs8);
const double clhs86 = K_darcy*N[2];
const double clhs87 = N[2]*rho;
const double clhs88 = bdf0*clhs87;
const double clhs89 = clhs87*volume_error_ratio;
const double clhs90 = clhs85 + clhs86 + clhs88 + clhs89;
const double clhs91 = clhs24*clhs90;
const double clhs92 = N[0]*clhs85 + clhs21*clhs90 - clhs23*clhs90 + clhs27*clhs91;
const double clhs93 = C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + clhs78;
const double clhs94 = C(1,3)*DN(2,1);
const double clhs95 = C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + clhs94;
const double clhs96 = C(3,5)*DN(2,0);
const double clhs97 = C(4,5)*DN(2,2);
const double clhs98 = C(1,5)*DN(2,1) + clhs96 + clhs97;
const double clhs99 = DN(2,1)*clhs35;
const double clhs100 = C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + clhs80;
const double clhs101 = C(3,4)*DN(2,1);
const double clhs102 = C(2,3)*DN(2,2) + clhs101 + clhs96;
const double clhs103 = C(2,5)*DN(2,2);
const double clhs104 = C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + clhs103;
const double clhs105 = DN(2,2)*clhs35;
const double clhs106 = DN(0,0)*N[2];
const double clhs107 = DN(2,0)*N[0];
const double clhs108 = C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double clhs109 = C(0,3)*DN(3,0);
const double clhs110 = C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + clhs109;
const double clhs111 = C(0,5)*DN(3,0);
const double clhs112 = C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + clhs111;
const double clhs113 = DN(0,0)*DN(3,0);
const double clhs114 = N[3]*clhs14 + N[3]*clhs16 + N[3]*clhs17;
const double clhs115 = clhs10*clhs113 + clhs114;
const double clhs116 = rho*(DN(3,0)*clhs6 + DN(3,1)*clhs7 + DN(3,2)*clhs8);
const double clhs117 = K_darcy*N[3];
const double clhs118 = N[3]*rho;
const double clhs119 = bdf0*clhs118;
const double clhs120 = clhs118*volume_error_ratio;
const double clhs121 = clhs116 + clhs117 + clhs119 + clhs120;
const double clhs122 = clhs121*clhs24;
const double clhs123 = N[0]*clhs116 + clhs121*clhs21 - clhs121*clhs23 + clhs122*clhs27;
const double clhs124 = C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + clhs109;
const double clhs125 = C(1,3)*DN(3,1);
const double clhs126 = C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + clhs125;
const double clhs127 = C(3,5)*DN(3,0);
const double clhs128 = C(4,5)*DN(3,2);
const double clhs129 = C(1,5)*DN(3,1) + clhs127 + clhs128;
const double clhs130 = DN(3,1)*clhs35;
const double clhs131 = C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + clhs111;
const double clhs132 = C(3,4)*DN(3,1);
const double clhs133 = C(2,3)*DN(3,2) + clhs127 + clhs132;
const double clhs134 = C(2,5)*DN(3,2);
const double clhs135 = C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + clhs134;
const double clhs136 = DN(3,2)*clhs35;
const double clhs137 = DN(0,0)*N[3];
const double clhs138 = DN(3,0)*N[0];
const double clhs139 = C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + clhs30;
const double clhs140 = C(0,4)*DN(0,0) + clhs33 + clhs38;
const double clhs141 = C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double clhs142 = C(1,4)*DN(0,1);
const double clhs143 = C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + clhs142;
const double clhs144 = pow(DN(0,1), 2);
const double clhs145 = C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + clhs142;
const double clhs146 = C(2,4)*DN(0,2);
const double clhs147 = C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + clhs146;
const double clhs148 = DN(0,1)*clhs10;
const double clhs149 = DN(0,2)*clhs148;
const double clhs150 = C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + clhs62;
const double clhs151 = C(0,4)*DN(1,0) + clhs65 + clhs69;
const double clhs152 = DN(1,0)*clhs148;
const double clhs153 = C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double clhs154 = C(1,4)*DN(1,1);
const double clhs155 = C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + clhs154;
const double clhs156 = DN(0,1)*DN(1,1);
const double clhs157 = clhs10*clhs156;
const double clhs158 = clhs51 + clhs60;
const double clhs159 = C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + clhs154;
const double clhs160 = C(2,4)*DN(1,2);
const double clhs161 = C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + clhs160;
const double clhs162 = DN(1,2)*clhs148;
const double clhs163 = DN(0,1)*N[1];
const double clhs164 = DN(1,1)*N[0];
const double clhs165 = C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + clhs94;
const double clhs166 = C(0,4)*DN(2,0) + clhs101 + clhs97;
const double clhs167 = DN(2,0)*clhs148;
const double clhs168 = C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double clhs169 = C(1,4)*DN(2,1);
const double clhs170 = C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + clhs169;
const double clhs171 = DN(0,1)*DN(2,1);
const double clhs172 = clhs10*clhs171;
const double clhs173 = clhs83 + clhs92;
const double clhs174 = C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + clhs169;
const double clhs175 = C(2,4)*DN(2,2);
const double clhs176 = C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + clhs175;
const double clhs177 = DN(2,2)*clhs148;
const double clhs178 = DN(0,1)*N[2];
const double clhs179 = DN(2,1)*N[0];
const double clhs180 = C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + clhs125;
const double clhs181 = C(0,4)*DN(3,0) + clhs128 + clhs132;
const double clhs182 = DN(3,0)*clhs148;
const double clhs183 = C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double clhs184 = C(1,4)*DN(3,1);
const double clhs185 = C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + clhs184;
const double clhs186 = DN(0,1)*DN(3,1);
const double clhs187 = clhs10*clhs186;
const double clhs188 = clhs114 + clhs123;
const double clhs189 = C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + clhs184;
const double clhs190 = C(2,4)*DN(3,2);
const double clhs191 = C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + clhs190;
const double clhs192 = DN(3,2)*clhs148;
const double clhs193 = DN(0,1)*N[3];
const double clhs194 = DN(3,1)*N[0];
const double clhs195 = C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + clhs40;
const double clhs196 = C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + clhs146;
const double clhs197 = C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double clhs198 = pow(DN(0,2), 2);
const double clhs199 = C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + clhs71;
const double clhs200 = DN(0,2)*clhs10;
const double clhs201 = DN(1,0)*clhs200;
const double clhs202 = C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + clhs160;
const double clhs203 = DN(1,1)*clhs200;
const double clhs204 = C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double clhs205 = DN(0,2)*DN(1,2);
const double clhs206 = clhs10*clhs205;
const double clhs207 = DN(0,2)*N[1];
const double clhs208 = DN(1,2)*N[0];
const double clhs209 = C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + clhs103;
const double clhs210 = DN(2,0)*clhs200;
const double clhs211 = C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + clhs175;
const double clhs212 = DN(2,1)*clhs200;
const double clhs213 = C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double clhs214 = DN(0,2)*DN(2,2);
const double clhs215 = clhs10*clhs214;
const double clhs216 = DN(0,2)*N[2];
const double clhs217 = DN(2,2)*N[0];
const double clhs218 = C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + clhs134;
const double clhs219 = DN(3,0)*clhs200;
const double clhs220 = C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + clhs190;
const double clhs221 = DN(3,1)*clhs200;
const double clhs222 = C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double clhs223 = DN(0,2)*DN(3,2);
const double clhs224 = clhs10*clhs223;
const double clhs225 = DN(0,2)*N[3];
const double clhs226 = DN(3,2)*N[0];
const double clhs227 = N[0] + clhs19*(1.0*clhs16 + 1.0*clhs17 + clhs20 + clhs22);
const double clhs228 = clhs24*(clhs156 + clhs205 + clhs50);
const double clhs229 = clhs24*(clhs171 + clhs214 + clhs82);
const double clhs230 = clhs24*(clhs113 + clhs186 + clhs223);
const double clhs231 = clhs24*clhs53;
const double clhs232 = clhs24*clhs54;
const double clhs233 = clhs26*clhs55;
const double clhs234 = N[1]*clhs12 + clhs18*clhs231 - clhs18*clhs232 + clhs233*clhs25;
const double clhs235 = pow(DN(1,0), 2);
const double clhs236 = pow(N[1], 2);
const double clhs237 = clhs236*rho;
const double clhs238 = K_darcy*clhs236 + N[1]*clhs53 + bdf0*clhs237 + clhs233*clhs59 + clhs237*volume_error_ratio + clhs53*clhs59 - clhs54*clhs59;
const double clhs239 = DN(1,0)*clhs10;
const double clhs240 = DN(1,1)*clhs239;
const double clhs241 = DN(1,2)*clhs239;
const double clhs242 = -N[1] + clhs231 - clhs232 + clhs43*clhs55;
const double clhs243 = DN(1,0)*DN(2,0);
const double clhs244 = N[2]*clhs54 + N[2]*clhs56 + N[2]*clhs57;
const double clhs245 = clhs10*clhs243 + clhs244;
const double clhs246 = N[1]*clhs85 + clhs233*clhs91 + clhs53*clhs91 - clhs54*clhs91;
const double clhs247 = DN(2,1)*clhs239;
const double clhs248 = DN(2,2)*clhs239;
const double clhs249 = DN(1,0)*N[2];
const double clhs250 = DN(2,0)*N[1];
const double clhs251 = DN(1,0)*DN(3,0);
const double clhs252 = N[3]*clhs54 + N[3]*clhs56 + N[3]*clhs57;
const double clhs253 = clhs10*clhs251 + clhs252;
const double clhs254 = N[1]*clhs116 + clhs122*clhs233 + clhs122*clhs53 - clhs122*clhs54;
const double clhs255 = DN(3,1)*clhs239;
const double clhs256 = DN(3,2)*clhs239;
const double clhs257 = DN(1,0)*N[3];
const double clhs258 = DN(3,0)*N[1];
const double clhs259 = clhs234 + clhs51;
const double clhs260 = pow(DN(1,1), 2);
const double clhs261 = DN(1,1)*clhs10;
const double clhs262 = DN(1,2)*clhs261;
const double clhs263 = DN(2,0)*clhs261;
const double clhs264 = DN(1,1)*DN(2,1);
const double clhs265 = clhs10*clhs264;
const double clhs266 = clhs244 + clhs246;
const double clhs267 = DN(2,2)*clhs261;
const double clhs268 = DN(1,1)*N[2];
const double clhs269 = DN(2,1)*N[1];
const double clhs270 = DN(3,0)*clhs261;
const double clhs271 = DN(1,1)*DN(3,1);
const double clhs272 = clhs10*clhs271;
const double clhs273 = clhs252 + clhs254;
const double clhs274 = DN(3,2)*clhs261;
const double clhs275 = DN(1,1)*N[3];
const double clhs276 = DN(3,1)*N[1];
const double clhs277 = pow(DN(1,2), 2);
const double clhs278 = DN(1,2)*clhs10;
const double clhs279 = DN(2,0)*clhs278;
const double clhs280 = DN(2,1)*clhs278;
const double clhs281 = DN(1,2)*DN(2,2);
const double clhs282 = clhs10*clhs281;
const double clhs283 = DN(1,2)*N[2];
const double clhs284 = DN(2,2)*N[1];
const double clhs285 = DN(3,0)*clhs278;
const double clhs286 = DN(3,1)*clhs278;
const double clhs287 = DN(1,2)*DN(3,2);
const double clhs288 = clhs10*clhs287;
const double clhs289 = DN(1,2)*N[3];
const double clhs290 = DN(3,2)*N[1];
const double clhs291 = N[1] + clhs19*(1.0*clhs53 + 1.0*clhs54 + 1.0*clhs56 + 1.0*clhs57);
const double clhs292 = clhs24*(clhs243 + clhs264 + clhs281);
const double clhs293 = clhs24*(clhs251 + clhs271 + clhs287);
const double clhs294 = clhs26*clhs87;
const double clhs295 = N[2]*clhs12 + clhs25*clhs294 + clhs25*clhs85 - clhs25*clhs86;
const double clhs296 = clhs24*clhs86;
const double clhs297 = clhs24*clhs85;
const double clhs298 = N[2]*clhs53 + clhs294*clhs59 + clhs59*clhs85 - clhs59*clhs86;
const double clhs299 = pow(DN(2,0), 2);
const double clhs300 = pow(N[2], 2);
const double clhs301 = clhs300*rho;
const double clhs302 = K_darcy*clhs300 + N[2]*clhs85 + bdf0*clhs301 + clhs294*clhs91 + clhs301*volume_error_ratio + clhs85*clhs91 - clhs86*clhs91;
const double clhs303 = DN(2,0)*clhs10;
const double clhs304 = DN(2,1)*clhs303;
const double clhs305 = DN(2,2)*clhs303;
const double clhs306 = -N[2] - clhs296 + clhs297 + clhs43*clhs87;
const double clhs307 = DN(2,0)*DN(3,0);
const double clhs308 = N[3]*clhs86 + N[3]*clhs88 + N[3]*clhs89;
const double clhs309 = clhs10*clhs307 + clhs308;
const double clhs310 = N[2]*clhs116 + clhs122*clhs294 + clhs122*clhs85 - clhs122*clhs86;
const double clhs311 = DN(3,1)*clhs303;
const double clhs312 = DN(3,2)*clhs303;
const double clhs313 = DN(2,0)*N[3];
const double clhs314 = DN(3,0)*N[2];
const double clhs315 = clhs295 + clhs83;
const double clhs316 = clhs244 + clhs298;
const double clhs317 = pow(DN(2,1), 2);
const double clhs318 = DN(2,1)*clhs10;
const double clhs319 = DN(2,2)*clhs318;
const double clhs320 = DN(3,0)*clhs318;
const double clhs321 = DN(2,1)*DN(3,1);
const double clhs322 = clhs10*clhs321;
const double clhs323 = clhs308 + clhs310;
const double clhs324 = DN(3,2)*clhs318;
const double clhs325 = DN(2,1)*N[3];
const double clhs326 = DN(3,1)*N[2];
const double clhs327 = pow(DN(2,2), 2);
const double clhs328 = DN(2,2)*clhs10;
const double clhs329 = DN(3,0)*clhs328;
const double clhs330 = DN(3,1)*clhs328;
const double clhs331 = DN(2,2)*DN(3,2);
const double clhs332 = clhs10*clhs331;
const double clhs333 = DN(2,2)*N[3];
const double clhs334 = DN(3,2)*N[2];
const double clhs335 = N[2] + clhs19*(1.0*clhs85 + 1.0*clhs86 + 1.0*clhs88 + 1.0*clhs89);
const double clhs336 = clhs24*(clhs307 + clhs321 + clhs331);
const double clhs337 = clhs118*clhs26;
const double clhs338 = N[3]*clhs12 + clhs116*clhs25 - clhs117*clhs25 + clhs25*clhs337;
const double clhs339 = clhs117*clhs24;
const double clhs340 = clhs116*clhs24;
const double clhs341 = N[3]*clhs53 + clhs116*clhs59 - clhs117*clhs59 + clhs337*clhs59;
const double clhs342 = N[3]*clhs85 + clhs116*clhs91 - clhs117*clhs91 + clhs337*clhs91;
const double clhs343 = pow(DN(3,0), 2);
const double clhs344 = pow(N[3], 2);
const double clhs345 = clhs344*rho;
const double clhs346 = K_darcy*clhs344 + N[3]*clhs116 + bdf0*clhs345 + clhs116*clhs122 - clhs117*clhs122 + clhs122*clhs337 + clhs345*volume_error_ratio;
const double clhs347 = DN(3,0)*clhs10;
const double clhs348 = DN(3,1)*clhs347;
const double clhs349 = DN(3,2)*clhs347;
const double clhs350 = -N[3] + clhs118*clhs43 - clhs339 + clhs340;
const double clhs351 = clhs114 + clhs338;
const double clhs352 = clhs252 + clhs341;
const double clhs353 = clhs308 + clhs342;
const double clhs354 = pow(DN(3,1), 2);
const double clhs355 = DN(3,1)*DN(3,2)*clhs10;
const double clhs356 = pow(DN(3,2), 2);
const double clhs357 = N[3] + clhs19*(1.0*clhs116 + 1.0*clhs117 + 1.0*clhs119 + 1.0*clhs120);
lhs(0,0)=DN(0,0)*clhs0 + DN(0,1)*clhs2 + DN(0,2)*clhs4 + clhs10*clhs5 + clhs28;
lhs(0,1)=DN(0,0)*clhs29 + DN(0,1)*clhs31 + DN(0,2)*clhs34 + clhs36;
lhs(0,2)=DN(0,0)*clhs37 + DN(0,1)*clhs39 + DN(0,2)*clhs41 + clhs42;
lhs(0,3)=DN(0,0)*clhs44;
lhs(0,4)=DN(0,0)*clhs45 + DN(0,1)*clhs47 + DN(0,2)*clhs49 + clhs52 + clhs60;
lhs(0,5)=DN(0,0)*clhs61 + DN(0,1)*clhs63 + DN(0,2)*clhs66 + clhs67;
lhs(0,6)=DN(0,0)*clhs68 + DN(0,1)*clhs70 + DN(0,2)*clhs72 + clhs73;
lhs(0,7)=DN(1,0)*clhs21 - DN(1,0)*clhs23 - clhs74 + clhs75*clhs76;
lhs(0,8)=DN(0,0)*clhs77 + DN(0,1)*clhs79 + DN(0,2)*clhs81 + clhs84 + clhs92;
lhs(0,9)=DN(0,0)*clhs93 + DN(0,1)*clhs95 + DN(0,2)*clhs98 + clhs99;
lhs(0,10)=DN(0,0)*clhs100 + DN(0,1)*clhs102 + DN(0,2)*clhs104 + clhs105;
lhs(0,11)=DN(2,0)*clhs21 - DN(2,0)*clhs23 - clhs106 + clhs107*clhs76;
lhs(0,12)=DN(0,0)*clhs108 + DN(0,1)*clhs110 + DN(0,2)*clhs112 + clhs115 + clhs123;
lhs(0,13)=DN(0,0)*clhs124 + DN(0,1)*clhs126 + DN(0,2)*clhs129 + clhs130;
lhs(0,14)=DN(0,0)*clhs131 + DN(0,1)*clhs133 + DN(0,2)*clhs135 + clhs136;
lhs(0,15)=DN(3,0)*clhs21 - DN(3,0)*clhs23 - clhs137 + clhs138*clhs76;
lhs(1,0)=DN(0,0)*clhs2 + DN(0,1)*clhs139 + DN(0,2)*clhs140 + clhs36;
lhs(1,1)=DN(0,0)*clhs31 + DN(0,1)*clhs141 + DN(0,2)*clhs143 + clhs10*clhs144 + clhs28;
lhs(1,2)=DN(0,0)*clhs39 + DN(0,1)*clhs145 + DN(0,2)*clhs147 + clhs149;
lhs(1,3)=DN(0,1)*clhs44;
lhs(1,4)=DN(0,0)*clhs47 + DN(0,1)*clhs150 + DN(0,2)*clhs151 + clhs152;
lhs(1,5)=DN(0,0)*clhs63 + DN(0,1)*clhs153 + DN(0,2)*clhs155 + clhs157 + clhs158;
lhs(1,6)=DN(0,0)*clhs70 + DN(0,1)*clhs159 + DN(0,2)*clhs161 + clhs162;
lhs(1,7)=DN(1,1)*clhs21 - DN(1,1)*clhs23 - clhs163 + clhs164*clhs76;
lhs(1,8)=DN(0,0)*clhs79 + DN(0,1)*clhs165 + DN(0,2)*clhs166 + clhs167;
lhs(1,9)=DN(0,0)*clhs95 + DN(0,1)*clhs168 + DN(0,2)*clhs170 + clhs172 + clhs173;
lhs(1,10)=DN(0,0)*clhs102 + DN(0,1)*clhs174 + DN(0,2)*clhs176 + clhs177;
lhs(1,11)=DN(2,1)*clhs21 - DN(2,1)*clhs23 - clhs178 + clhs179*clhs76;
lhs(1,12)=DN(0,0)*clhs110 + DN(0,1)*clhs180 + DN(0,2)*clhs181 + clhs182;
lhs(1,13)=DN(0,0)*clhs126 + DN(0,1)*clhs183 + DN(0,2)*clhs185 + clhs187 + clhs188;
lhs(1,14)=DN(0,0)*clhs133 + DN(0,1)*clhs189 + DN(0,2)*clhs191 + clhs192;
lhs(1,15)=DN(3,1)*clhs21 - DN(3,1)*clhs23 - clhs193 + clhs194*clhs76;
lhs(2,0)=DN(0,0)*clhs4 + DN(0,1)*clhs140 + DN(0,2)*clhs195 + clhs42;
lhs(2,1)=DN(0,0)*clhs34 + DN(0,1)*clhs143 + DN(0,2)*clhs196 + clhs149;
lhs(2,2)=DN(0,0)*clhs41 + DN(0,1)*clhs147 + DN(0,2)*clhs197 + clhs10*clhs198 + clhs28;
lhs(2,3)=DN(0,2)*clhs44;
lhs(2,4)=DN(0,0)*clhs49 + DN(0,1)*clhs151 + DN(0,2)*clhs199 + clhs201;
lhs(2,5)=DN(0,0)*clhs66 + DN(0,1)*clhs155 + DN(0,2)*clhs202 + clhs203;
lhs(2,6)=DN(0,0)*clhs72 + DN(0,1)*clhs161 + DN(0,2)*clhs204 + clhs158 + clhs206;
lhs(2,7)=DN(1,2)*clhs21 - DN(1,2)*clhs23 - clhs207 + clhs208*clhs76;
lhs(2,8)=DN(0,0)*clhs81 + DN(0,1)*clhs166 + DN(0,2)*clhs209 + clhs210;
lhs(2,9)=DN(0,0)*clhs98 + DN(0,1)*clhs170 + DN(0,2)*clhs211 + clhs212;
lhs(2,10)=DN(0,0)*clhs104 + DN(0,1)*clhs176 + DN(0,2)*clhs213 + clhs173 + clhs215;
lhs(2,11)=DN(2,2)*clhs21 - DN(2,2)*clhs23 - clhs216 + clhs217*clhs76;
lhs(2,12)=DN(0,0)*clhs112 + DN(0,1)*clhs181 + DN(0,2)*clhs218 + clhs219;
lhs(2,13)=DN(0,0)*clhs129 + DN(0,1)*clhs185 + DN(0,2)*clhs220 + clhs221;
lhs(2,14)=DN(0,0)*clhs135 + DN(0,1)*clhs191 + DN(0,2)*clhs222 + clhs188 + clhs224;
lhs(2,15)=DN(3,2)*clhs21 - DN(3,2)*clhs23 - clhs225 + clhs226*clhs76;
lhs(3,0)=DN(0,0)*clhs227;
lhs(3,1)=DN(0,1)*clhs227;
lhs(3,2)=DN(0,2)*clhs227;
lhs(3,3)=clhs24*(clhs144 + clhs198 + clhs5);
lhs(3,4)=DN(0,0)*clhs59 + clhs75;
lhs(3,5)=DN(0,1)*clhs59 + clhs164;
lhs(3,6)=DN(0,2)*clhs59 + clhs208;
lhs(3,7)=clhs228;
lhs(3,8)=DN(0,0)*clhs91 + clhs107;
lhs(3,9)=DN(0,1)*clhs91 + clhs179;
lhs(3,10)=DN(0,2)*clhs91 + clhs217;
lhs(3,11)=clhs229;
lhs(3,12)=DN(0,0)*clhs122 + clhs138;
lhs(3,13)=DN(0,1)*clhs122 + clhs194;
lhs(3,14)=DN(0,2)*clhs122 + clhs226;
lhs(3,15)=clhs230;
lhs(4,0)=DN(1,0)*clhs0 + DN(1,1)*clhs2 + DN(1,2)*clhs4 + clhs234 + clhs52;
lhs(4,1)=DN(1,0)*clhs29 + DN(1,1)*clhs31 + DN(1,2)*clhs34 + clhs152;
lhs(4,2)=DN(1,0)*clhs37 + DN(1,1)*clhs39 + DN(1,2)*clhs41 + clhs201;
lhs(4,3)=DN(0,0)*clhs231 - DN(0,0)*clhs232 + clhs74*clhs76 - clhs75;
lhs(4,4)=DN(1,0)*clhs45 + DN(1,1)*clhs47 + DN(1,2)*clhs49 + clhs10*clhs235 + clhs238;
lhs(4,5)=DN(1,0)*clhs61 + DN(1,1)*clhs63 + DN(1,2)*clhs66 + clhs240;
lhs(4,6)=DN(1,0)*clhs68 + DN(1,1)*clhs70 + DN(1,2)*clhs72 + clhs241;
lhs(4,7)=DN(1,0)*clhs242;
lhs(4,8)=DN(1,0)*clhs77 + DN(1,1)*clhs79 + DN(1,2)*clhs81 + clhs245 + clhs246;
lhs(4,9)=DN(1,0)*clhs93 + DN(1,1)*clhs95 + DN(1,2)*clhs98 + clhs247;
lhs(4,10)=DN(1,0)*clhs100 + DN(1,1)*clhs102 + DN(1,2)*clhs104 + clhs248;
lhs(4,11)=DN(2,0)*clhs231 - DN(2,0)*clhs232 - clhs249 + clhs250*clhs76;
lhs(4,12)=DN(1,0)*clhs108 + DN(1,1)*clhs110 + DN(1,2)*clhs112 + clhs253 + clhs254;
lhs(4,13)=DN(1,0)*clhs124 + DN(1,1)*clhs126 + DN(1,2)*clhs129 + clhs255;
lhs(4,14)=DN(1,0)*clhs131 + DN(1,1)*clhs133 + DN(1,2)*clhs135 + clhs256;
lhs(4,15)=DN(3,0)*clhs231 - DN(3,0)*clhs232 - clhs257 + clhs258*clhs76;
lhs(5,0)=DN(1,0)*clhs2 + DN(1,1)*clhs139 + DN(1,2)*clhs140 + clhs67;
lhs(5,1)=DN(1,0)*clhs31 + DN(1,1)*clhs141 + DN(1,2)*clhs143 + clhs157 + clhs259;
lhs(5,2)=DN(1,0)*clhs39 + DN(1,1)*clhs145 + DN(1,2)*clhs147 + clhs203;
lhs(5,3)=DN(0,1)*clhs231 - DN(0,1)*clhs232 + clhs163*clhs76 - clhs164;
lhs(5,4)=DN(1,0)*clhs47 + DN(1,1)*clhs150 + DN(1,2)*clhs151 + clhs240;
lhs(5,5)=DN(1,0)*clhs63 + DN(1,1)*clhs153 + DN(1,2)*clhs155 + clhs10*clhs260 + clhs238;
lhs(5,6)=DN(1,0)*clhs70 + DN(1,1)*clhs159 + DN(1,2)*clhs161 + clhs262;
lhs(5,7)=DN(1,1)*clhs242;
lhs(5,8)=DN(1,0)*clhs79 + DN(1,1)*clhs165 + DN(1,2)*clhs166 + clhs263;
lhs(5,9)=DN(1,0)*clhs95 + DN(1,1)*clhs168 + DN(1,2)*clhs170 + clhs265 + clhs266;
lhs(5,10)=DN(1,0)*clhs102 + DN(1,1)*clhs174 + DN(1,2)*clhs176 + clhs267;
lhs(5,11)=DN(2,1)*clhs231 - DN(2,1)*clhs232 - clhs268 + clhs269*clhs76;
lhs(5,12)=DN(1,0)*clhs110 + DN(1,1)*clhs180 + DN(1,2)*clhs181 + clhs270;
lhs(5,13)=DN(1,0)*clhs126 + DN(1,1)*clhs183 + DN(1,2)*clhs185 + clhs272 + clhs273;
lhs(5,14)=DN(1,0)*clhs133 + DN(1,1)*clhs189 + DN(1,2)*clhs191 + clhs274;
lhs(5,15)=DN(3,1)*clhs231 - DN(3,1)*clhs232 - clhs275 + clhs276*clhs76;
lhs(6,0)=DN(1,0)*clhs4 + DN(1,1)*clhs140 + DN(1,2)*clhs195 + clhs73;
lhs(6,1)=DN(1,0)*clhs34 + DN(1,1)*clhs143 + DN(1,2)*clhs196 + clhs162;
lhs(6,2)=DN(1,0)*clhs41 + DN(1,1)*clhs147 + DN(1,2)*clhs197 + clhs206 + clhs259;
lhs(6,3)=DN(0,2)*clhs231 - DN(0,2)*clhs232 + clhs207*clhs76 - clhs208;
lhs(6,4)=DN(1,0)*clhs49 + DN(1,1)*clhs151 + DN(1,2)*clhs199 + clhs241;
lhs(6,5)=DN(1,0)*clhs66 + DN(1,1)*clhs155 + DN(1,2)*clhs202 + clhs262;
lhs(6,6)=DN(1,0)*clhs72 + DN(1,1)*clhs161 + DN(1,2)*clhs204 + clhs10*clhs277 + clhs238;
lhs(6,7)=DN(1,2)*clhs242;
lhs(6,8)=DN(1,0)*clhs81 + DN(1,1)*clhs166 + DN(1,2)*clhs209 + clhs279;
lhs(6,9)=DN(1,0)*clhs98 + DN(1,1)*clhs170 + DN(1,2)*clhs211 + clhs280;
lhs(6,10)=DN(1,0)*clhs104 + DN(1,1)*clhs176 + DN(1,2)*clhs213 + clhs266 + clhs282;
lhs(6,11)=DN(2,2)*clhs231 - DN(2,2)*clhs232 - clhs283 + clhs284*clhs76;
lhs(6,12)=DN(1,0)*clhs112 + DN(1,1)*clhs181 + DN(1,2)*clhs218 + clhs285;
lhs(6,13)=DN(1,0)*clhs129 + DN(1,1)*clhs185 + DN(1,2)*clhs220 + clhs286;
lhs(6,14)=DN(1,0)*clhs135 + DN(1,1)*clhs191 + DN(1,2)*clhs222 + clhs273 + clhs288;
lhs(6,15)=DN(3,2)*clhs231 - DN(3,2)*clhs232 - clhs289 + clhs290*clhs76;
lhs(7,0)=DN(1,0)*clhs25 + clhs74;
lhs(7,1)=DN(1,1)*clhs25 + clhs163;
lhs(7,2)=DN(1,2)*clhs25 + clhs207;
lhs(7,3)=clhs228;
lhs(7,4)=DN(1,0)*clhs291;
lhs(7,5)=DN(1,1)*clhs291;
lhs(7,6)=DN(1,2)*clhs291;
lhs(7,7)=clhs24*(clhs235 + clhs260 + clhs277);
lhs(7,8)=DN(1,0)*clhs91 + clhs250;
lhs(7,9)=DN(1,1)*clhs91 + clhs269;
lhs(7,10)=DN(1,2)*clhs91 + clhs284;
lhs(7,11)=clhs292;
lhs(7,12)=DN(1,0)*clhs122 + clhs258;
lhs(7,13)=DN(1,1)*clhs122 + clhs276;
lhs(7,14)=DN(1,2)*clhs122 + clhs290;
lhs(7,15)=clhs293;
lhs(8,0)=DN(2,0)*clhs0 + DN(2,1)*clhs2 + DN(2,2)*clhs4 + clhs295 + clhs84;
lhs(8,1)=DN(2,0)*clhs29 + DN(2,1)*clhs31 + DN(2,2)*clhs34 + clhs167;
lhs(8,2)=DN(2,0)*clhs37 + DN(2,1)*clhs39 + DN(2,2)*clhs41 + clhs210;
lhs(8,3)=-DN(0,0)*clhs296 + DN(0,0)*clhs297 + clhs106*clhs76 - clhs107;
lhs(8,4)=DN(2,0)*clhs45 + DN(2,1)*clhs47 + DN(2,2)*clhs49 + clhs245 + clhs298;
lhs(8,5)=DN(2,0)*clhs61 + DN(2,1)*clhs63 + DN(2,2)*clhs66 + clhs263;
lhs(8,6)=DN(2,0)*clhs68 + DN(2,1)*clhs70 + DN(2,2)*clhs72 + clhs279;
lhs(8,7)=-DN(1,0)*clhs296 + DN(1,0)*clhs297 + clhs249*clhs76 - clhs250;
lhs(8,8)=DN(2,0)*clhs77 + DN(2,1)*clhs79 + DN(2,2)*clhs81 + clhs10*clhs299 + clhs302;
lhs(8,9)=DN(2,0)*clhs93 + DN(2,1)*clhs95 + DN(2,2)*clhs98 + clhs304;
lhs(8,10)=DN(2,0)*clhs100 + DN(2,1)*clhs102 + DN(2,2)*clhs104 + clhs305;
lhs(8,11)=DN(2,0)*clhs306;
lhs(8,12)=DN(2,0)*clhs108 + DN(2,1)*clhs110 + DN(2,2)*clhs112 + clhs309 + clhs310;
lhs(8,13)=DN(2,0)*clhs124 + DN(2,1)*clhs126 + DN(2,2)*clhs129 + clhs311;
lhs(8,14)=DN(2,0)*clhs131 + DN(2,1)*clhs133 + DN(2,2)*clhs135 + clhs312;
lhs(8,15)=-DN(3,0)*clhs296 + DN(3,0)*clhs297 - clhs313 + clhs314*clhs76;
lhs(9,0)=DN(2,0)*clhs2 + DN(2,1)*clhs139 + DN(2,2)*clhs140 + clhs99;
lhs(9,1)=DN(2,0)*clhs31 + DN(2,1)*clhs141 + DN(2,2)*clhs143 + clhs172 + clhs315;
lhs(9,2)=DN(2,0)*clhs39 + DN(2,1)*clhs145 + DN(2,2)*clhs147 + clhs212;
lhs(9,3)=-DN(0,1)*clhs296 + DN(0,1)*clhs297 + clhs178*clhs76 - clhs179;
lhs(9,4)=DN(2,0)*clhs47 + DN(2,1)*clhs150 + DN(2,2)*clhs151 + clhs247;
lhs(9,5)=DN(2,0)*clhs63 + DN(2,1)*clhs153 + DN(2,2)*clhs155 + clhs265 + clhs316;
lhs(9,6)=DN(2,0)*clhs70 + DN(2,1)*clhs159 + DN(2,2)*clhs161 + clhs280;
lhs(9,7)=-DN(1,1)*clhs296 + DN(1,1)*clhs297 + clhs268*clhs76 - clhs269;
lhs(9,8)=DN(2,0)*clhs79 + DN(2,1)*clhs165 + DN(2,2)*clhs166 + clhs304;
lhs(9,9)=DN(2,0)*clhs95 + DN(2,1)*clhs168 + DN(2,2)*clhs170 + clhs10*clhs317 + clhs302;
lhs(9,10)=DN(2,0)*clhs102 + DN(2,1)*clhs174 + DN(2,2)*clhs176 + clhs319;
lhs(9,11)=DN(2,1)*clhs306;
lhs(9,12)=DN(2,0)*clhs110 + DN(2,1)*clhs180 + DN(2,2)*clhs181 + clhs320;
lhs(9,13)=DN(2,0)*clhs126 + DN(2,1)*clhs183 + DN(2,2)*clhs185 + clhs322 + clhs323;
lhs(9,14)=DN(2,0)*clhs133 + DN(2,1)*clhs189 + DN(2,2)*clhs191 + clhs324;
lhs(9,15)=-DN(3,1)*clhs296 + DN(3,1)*clhs297 - clhs325 + clhs326*clhs76;
lhs(10,0)=DN(2,0)*clhs4 + DN(2,1)*clhs140 + DN(2,2)*clhs195 + clhs105;
lhs(10,1)=DN(2,0)*clhs34 + DN(2,1)*clhs143 + DN(2,2)*clhs196 + clhs177;
lhs(10,2)=DN(2,0)*clhs41 + DN(2,1)*clhs147 + DN(2,2)*clhs197 + clhs215 + clhs315;
lhs(10,3)=-DN(0,2)*clhs296 + DN(0,2)*clhs297 + clhs216*clhs76 - clhs217;
lhs(10,4)=DN(2,0)*clhs49 + DN(2,1)*clhs151 + DN(2,2)*clhs199 + clhs248;
lhs(10,5)=DN(2,0)*clhs66 + DN(2,1)*clhs155 + DN(2,2)*clhs202 + clhs267;
lhs(10,6)=DN(2,0)*clhs72 + DN(2,1)*clhs161 + DN(2,2)*clhs204 + clhs282 + clhs316;
lhs(10,7)=-DN(1,2)*clhs296 + DN(1,2)*clhs297 + clhs283*clhs76 - clhs284;
lhs(10,8)=DN(2,0)*clhs81 + DN(2,1)*clhs166 + DN(2,2)*clhs209 + clhs305;
lhs(10,9)=DN(2,0)*clhs98 + DN(2,1)*clhs170 + DN(2,2)*clhs211 + clhs319;
lhs(10,10)=DN(2,0)*clhs104 + DN(2,1)*clhs176 + DN(2,2)*clhs213 + clhs10*clhs327 + clhs302;
lhs(10,11)=DN(2,2)*clhs306;
lhs(10,12)=DN(2,0)*clhs112 + DN(2,1)*clhs181 + DN(2,2)*clhs218 + clhs329;
lhs(10,13)=DN(2,0)*clhs129 + DN(2,1)*clhs185 + DN(2,2)*clhs220 + clhs330;
lhs(10,14)=DN(2,0)*clhs135 + DN(2,1)*clhs191 + DN(2,2)*clhs222 + clhs323 + clhs332;
lhs(10,15)=-DN(3,2)*clhs296 + DN(3,2)*clhs297 - clhs333 + clhs334*clhs76;
lhs(11,0)=DN(2,0)*clhs25 + clhs106;
lhs(11,1)=DN(2,1)*clhs25 + clhs178;
lhs(11,2)=DN(2,2)*clhs25 + clhs216;
lhs(11,3)=clhs229;
lhs(11,4)=DN(2,0)*clhs59 + clhs249;
lhs(11,5)=DN(2,1)*clhs59 + clhs268;
lhs(11,6)=DN(2,2)*clhs59 + clhs283;
lhs(11,7)=clhs292;
lhs(11,8)=DN(2,0)*clhs335;
lhs(11,9)=DN(2,1)*clhs335;
lhs(11,10)=DN(2,2)*clhs335;
lhs(11,11)=clhs24*(clhs299 + clhs317 + clhs327);
lhs(11,12)=DN(2,0)*clhs122 + clhs314;
lhs(11,13)=DN(2,1)*clhs122 + clhs326;
lhs(11,14)=DN(2,2)*clhs122 + clhs334;
lhs(11,15)=clhs336;
lhs(12,0)=DN(3,0)*clhs0 + DN(3,1)*clhs2 + DN(3,2)*clhs4 + clhs115 + clhs338;
lhs(12,1)=DN(3,0)*clhs29 + DN(3,1)*clhs31 + DN(3,2)*clhs34 + clhs182;
lhs(12,2)=DN(3,0)*clhs37 + DN(3,1)*clhs39 + DN(3,2)*clhs41 + clhs219;
lhs(12,3)=-DN(0,0)*clhs339 + DN(0,0)*clhs340 + clhs137*clhs76 - clhs138;
lhs(12,4)=DN(3,0)*clhs45 + DN(3,1)*clhs47 + DN(3,2)*clhs49 + clhs253 + clhs341;
lhs(12,5)=DN(3,0)*clhs61 + DN(3,1)*clhs63 + DN(3,2)*clhs66 + clhs270;
lhs(12,6)=DN(3,0)*clhs68 + DN(3,1)*clhs70 + DN(3,2)*clhs72 + clhs285;
lhs(12,7)=-DN(1,0)*clhs339 + DN(1,0)*clhs340 + clhs257*clhs76 - clhs258;
lhs(12,8)=DN(3,0)*clhs77 + DN(3,1)*clhs79 + DN(3,2)*clhs81 + clhs309 + clhs342;
lhs(12,9)=DN(3,0)*clhs93 + DN(3,1)*clhs95 + DN(3,2)*clhs98 + clhs320;
lhs(12,10)=DN(3,0)*clhs100 + DN(3,1)*clhs102 + DN(3,2)*clhs104 + clhs329;
lhs(12,11)=-DN(2,0)*clhs339 + DN(2,0)*clhs340 + clhs313*clhs76 - clhs314;
lhs(12,12)=DN(3,0)*clhs108 + DN(3,1)*clhs110 + DN(3,2)*clhs112 + clhs10*clhs343 + clhs346;
lhs(12,13)=DN(3,0)*clhs124 + DN(3,1)*clhs126 + DN(3,2)*clhs129 + clhs348;
lhs(12,14)=DN(3,0)*clhs131 + DN(3,1)*clhs133 + DN(3,2)*clhs135 + clhs349;
lhs(12,15)=DN(3,0)*clhs350;
lhs(13,0)=DN(3,0)*clhs2 + DN(3,1)*clhs139 + DN(3,2)*clhs140 + clhs130;
lhs(13,1)=DN(3,0)*clhs31 + DN(3,1)*clhs141 + DN(3,2)*clhs143 + clhs187 + clhs351;
lhs(13,2)=DN(3,0)*clhs39 + DN(3,1)*clhs145 + DN(3,2)*clhs147 + clhs221;
lhs(13,3)=-DN(0,1)*clhs339 + DN(0,1)*clhs340 + clhs193*clhs76 - clhs194;
lhs(13,4)=DN(3,0)*clhs47 + DN(3,1)*clhs150 + DN(3,2)*clhs151 + clhs255;
lhs(13,5)=DN(3,0)*clhs63 + DN(3,1)*clhs153 + DN(3,2)*clhs155 + clhs272 + clhs352;
lhs(13,6)=DN(3,0)*clhs70 + DN(3,1)*clhs159 + DN(3,2)*clhs161 + clhs286;
lhs(13,7)=-DN(1,1)*clhs339 + DN(1,1)*clhs340 + clhs275*clhs76 - clhs276;
lhs(13,8)=DN(3,0)*clhs79 + DN(3,1)*clhs165 + DN(3,2)*clhs166 + clhs311;
lhs(13,9)=DN(3,0)*clhs95 + DN(3,1)*clhs168 + DN(3,2)*clhs170 + clhs322 + clhs353;
lhs(13,10)=DN(3,0)*clhs102 + DN(3,1)*clhs174 + DN(3,2)*clhs176 + clhs330;
lhs(13,11)=-DN(2,1)*clhs339 + DN(2,1)*clhs340 + clhs325*clhs76 - clhs326;
lhs(13,12)=DN(3,0)*clhs110 + DN(3,1)*clhs180 + DN(3,2)*clhs181 + clhs348;
lhs(13,13)=DN(3,0)*clhs126 + DN(3,1)*clhs183 + DN(3,2)*clhs185 + clhs10*clhs354 + clhs346;
lhs(13,14)=DN(3,0)*clhs133 + DN(3,1)*clhs189 + DN(3,2)*clhs191 + clhs355;
lhs(13,15)=DN(3,1)*clhs350;
lhs(14,0)=DN(3,0)*clhs4 + DN(3,1)*clhs140 + DN(3,2)*clhs195 + clhs136;
lhs(14,1)=DN(3,0)*clhs34 + DN(3,1)*clhs143 + DN(3,2)*clhs196 + clhs192;
lhs(14,2)=DN(3,0)*clhs41 + DN(3,1)*clhs147 + DN(3,2)*clhs197 + clhs224 + clhs351;
lhs(14,3)=-DN(0,2)*clhs339 + DN(0,2)*clhs340 + clhs225*clhs76 - clhs226;
lhs(14,4)=DN(3,0)*clhs49 + DN(3,1)*clhs151 + DN(3,2)*clhs199 + clhs256;
lhs(14,5)=DN(3,0)*clhs66 + DN(3,1)*clhs155 + DN(3,2)*clhs202 + clhs274;
lhs(14,6)=DN(3,0)*clhs72 + DN(3,1)*clhs161 + DN(3,2)*clhs204 + clhs288 + clhs352;
lhs(14,7)=-DN(1,2)*clhs339 + DN(1,2)*clhs340 + clhs289*clhs76 - clhs290;
lhs(14,8)=DN(3,0)*clhs81 + DN(3,1)*clhs166 + DN(3,2)*clhs209 + clhs312;
lhs(14,9)=DN(3,0)*clhs98 + DN(3,1)*clhs170 + DN(3,2)*clhs211 + clhs324;
lhs(14,10)=DN(3,0)*clhs104 + DN(3,1)*clhs176 + DN(3,2)*clhs213 + clhs332 + clhs353;
lhs(14,11)=-DN(2,2)*clhs339 + DN(2,2)*clhs340 + clhs333*clhs76 - clhs334;
lhs(14,12)=DN(3,0)*clhs112 + DN(3,1)*clhs181 + DN(3,2)*clhs218 + clhs349;
lhs(14,13)=DN(3,0)*clhs129 + DN(3,1)*clhs185 + DN(3,2)*clhs220 + clhs355;
lhs(14,14)=DN(3,0)*clhs135 + DN(3,1)*clhs191 + DN(3,2)*clhs222 + clhs10*clhs356 + clhs346;
lhs(14,15)=DN(3,2)*clhs350;
lhs(15,0)=DN(3,0)*clhs25 + clhs137;
lhs(15,1)=DN(3,1)*clhs25 + clhs193;
lhs(15,2)=DN(3,2)*clhs25 + clhs225;
lhs(15,3)=clhs230;
lhs(15,4)=DN(3,0)*clhs59 + clhs257;
lhs(15,5)=DN(3,1)*clhs59 + clhs275;
lhs(15,6)=DN(3,2)*clhs59 + clhs289;
lhs(15,7)=clhs293;
lhs(15,8)=DN(3,0)*clhs91 + clhs313;
lhs(15,9)=DN(3,1)*clhs91 + clhs325;
lhs(15,10)=DN(3,2)*clhs91 + clhs333;
lhs(15,11)=clhs336;
lhs(15,12)=DN(3,0)*clhs357;
lhs(15,13)=DN(3,1)*clhs357;
lhs(15,14)=DN(3,2)*clhs357;
lhs(15,15)=clhs24*(clhs343 + clhs354 + clhs356);


    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;
}

template <>
void TwoFluidNavierStokes<TwoFluidNavierStokesData<2, 3>>::ComputeGaussPointRHSContribution(
    TwoFluidNavierStokesData<2, 3> &rData,
    VectorType &rRHS)
{

    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double bdf1 = rData.bdf1;
    const double bdf2 = rData.bdf2;

    const double dyn_tau = rData.DynamicTau;
    const double K_darcy = rData.DarcyTerm;

    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vnn = rData.Velocity_OldStep2;
    const auto &vmesh = rData.MeshVelocity;
    const auto &vconv = v - vmesh;
    const auto &f = rData.BodyForce;
    const auto &p = rData.Pressure;
    const auto &stress = rData.ShearStress;

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    // Mass correction term

    const double volume_error_ratio = rData.VolumeError;

    auto &rhs = rData.rhs;

    const double crhs0 = N[0]*p[0] + N[1]*p[1] + N[2]*p[2];
const double crhs1 = rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0));
const double crhs2 = N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0);
const double crhs3 = K_darcy*crhs2;
const double crhs4 = rho*volume_error_ratio;
const double crhs5 = crhs2*crhs4;
const double crhs6 = rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)));
const double crhs7 = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
const double crhs8 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs9 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs10 = rho*(crhs7*crhs8 + crhs9*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0)));
const double crhs11 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
const double crhs12 = crhs11 + crhs7 - volume_error_ratio;
const double crhs13 = rho*stab_c2*sqrt(pow(crhs8, 2) + pow(crhs9, 2));
const double crhs14 = crhs12*(crhs13*h/stab_c1 + mu);
const double crhs15 = 1.0/(K_darcy + crhs13/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs16 = crhs15*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] - crhs1 + crhs10 + crhs3 + crhs5 + crhs6);
const double crhs17 = K_darcy*N[0];
const double crhs18 = rho*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1));
const double crhs19 = N[0]*crhs18;
const double crhs20 = rho*(DN(0,0)*crhs8 + DN(0,1)*crhs9);
const double crhs21 = rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1));
const double crhs22 = N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1);
const double crhs23 = K_darcy*crhs22;
const double crhs24 = crhs22*crhs4;
const double crhs25 = rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)));
const double crhs26 = rho*(crhs11*crhs9 + crhs8*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1)));
const double crhs27 = crhs15*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] - crhs21 + crhs23 + crhs24 + crhs25 + crhs26);
const double crhs28 = K_darcy*N[1];
const double crhs29 = N[1]*crhs18;
const double crhs30 = rho*(DN(1,0)*crhs8 + DN(1,1)*crhs9);
const double crhs31 = K_darcy*N[2];
const double crhs32 = N[2]*crhs18;
const double crhs33 = rho*(DN(2,0)*crhs8 + DN(2,1)*crhs9);
rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs14 - DN(0,0)*stress[0] - DN(0,1)*stress[2] + N[0]*crhs1 - N[0]*crhs10 - N[0]*crhs3 - N[0]*crhs5 - N[0]*crhs6 + crhs16*crhs17 - crhs16*crhs19 - crhs16*crhs20;
rhs[1]=-DN(0,0)*stress[2] + DN(0,1)*crhs0 - DN(0,1)*crhs14 - DN(0,1)*stress[1] + N[0]*crhs21 - N[0]*crhs23 - N[0]*crhs24 - N[0]*crhs25 - N[0]*crhs26 + crhs17*crhs27 - crhs19*crhs27 - crhs20*crhs27;
rhs[2]=-DN(0,0)*crhs16 - DN(0,1)*crhs27 - N[0]*crhs12;
rhs[3]=DN(1,0)*crhs0 - DN(1,0)*crhs14 - DN(1,0)*stress[0] - DN(1,1)*stress[2] + N[1]*crhs1 - N[1]*crhs10 - N[1]*crhs3 - N[1]*crhs5 - N[1]*crhs6 + crhs16*crhs28 - crhs16*crhs29 - crhs16*crhs30;
rhs[4]=-DN(1,0)*stress[2] + DN(1,1)*crhs0 - DN(1,1)*crhs14 - DN(1,1)*stress[1] + N[1]*crhs21 - N[1]*crhs23 - N[1]*crhs24 - N[1]*crhs25 - N[1]*crhs26 + crhs27*crhs28 - crhs27*crhs29 - crhs27*crhs30;
rhs[5]=-DN(1,0)*crhs16 - DN(1,1)*crhs27 - N[1]*crhs12;
rhs[6]=DN(2,0)*crhs0 - DN(2,0)*crhs14 - DN(2,0)*stress[0] - DN(2,1)*stress[2] + N[2]*crhs1 - N[2]*crhs10 - N[2]*crhs3 - N[2]*crhs5 - N[2]*crhs6 + crhs16*crhs31 - crhs16*crhs32 - crhs16*crhs33;
rhs[7]=-DN(2,0)*stress[2] + DN(2,1)*crhs0 - DN(2,1)*crhs14 - DN(2,1)*stress[1] + N[2]*crhs21 - N[2]*crhs23 - N[2]*crhs24 - N[2]*crhs25 - N[2]*crhs26 + crhs27*crhs31 - crhs27*crhs32 - crhs27*crhs33;
rhs[8]=-DN(2,0)*crhs16 - DN(2,1)*crhs27 - N[2]*crhs12;


    noalias(rRHS) += rData.Weight * rhs;
}

template <>
void TwoFluidNavierStokes<TwoFluidNavierStokesData<3, 4>>::ComputeGaussPointRHSContribution(
    TwoFluidNavierStokesData<3, 4> &rData,
    VectorType &rRHS)
{

    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double bdf1 = rData.bdf1;
    const double bdf2 = rData.bdf2;

    const double dyn_tau = rData.DynamicTau;
    const double K_darcy = rData.DarcyTerm;

    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vnn = rData.Velocity_OldStep2;
    const auto &vmesh = rData.MeshVelocity;
    const auto &vconv = v - vmesh;
    const auto &f = rData.BodyForce;
    const auto &p = rData.Pressure;
    const auto &stress = rData.ShearStress;

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    // Mass correction term
    const double volume_error_ratio = rData.VolumeError;

    auto &rhs = rData.rhs;

    const double crhs0 = N[0]*p[0] + N[1]*p[1] + N[2]*p[2] + N[3]*p[3];
const double crhs1 = rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0));
const double crhs2 = N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0);
const double crhs3 = K_darcy*crhs2;
const double crhs4 = rho*volume_error_ratio;
const double crhs5 = crhs2*crhs4;
const double crhs6 = rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0)));
const double crhs7 = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0);
const double crhs8 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs9 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs10 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs11 = rho*(crhs10*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0)) + crhs7*crhs8 + crhs9*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)));
const double crhs12 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
const double crhs13 = DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
const double crhs14 = crhs12 + crhs13 + crhs7 - volume_error_ratio;
const double crhs15 = rho*stab_c2*sqrt(pow(crhs10, 2) + pow(crhs8, 2) + pow(crhs9, 2));
const double crhs16 = crhs14*(crhs15*h/stab_c1 + mu);
const double crhs17 = 1.0/(K_darcy + crhs15/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs18 = crhs17*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] - crhs1 + crhs11 + crhs3 + crhs5 + crhs6);
const double crhs19 = K_darcy*N[0];
const double crhs20 = rho*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2));
const double crhs21 = N[0]*crhs20;
const double crhs22 = rho*(DN(0,0)*crhs8 + DN(0,1)*crhs9 + DN(0,2)*crhs10);
const double crhs23 = rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1));
const double crhs24 = N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1);
const double crhs25 = K_darcy*crhs24;
const double crhs26 = crhs24*crhs4;
const double crhs27 = rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1)));
const double crhs28 = rho*(crhs10*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1)) + crhs12*crhs9 + crhs8*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)));
const double crhs29 = crhs17*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] - crhs23 + crhs25 + crhs26 + crhs27 + crhs28);
const double crhs30 = rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2));
const double crhs31 = N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2);
const double crhs32 = K_darcy*crhs31;
const double crhs33 = crhs31*crhs4;
const double crhs34 = rho*(N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2)));
const double crhs35 = rho*(crhs10*crhs13 + crhs8*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + crhs9*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2)));
const double crhs36 = crhs17*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] - crhs30 + crhs32 + crhs33 + crhs34 + crhs35);
const double crhs37 = K_darcy*N[1];
const double crhs38 = N[1]*crhs20;
const double crhs39 = rho*(DN(1,0)*crhs8 + DN(1,1)*crhs9 + DN(1,2)*crhs10);
const double crhs40 = K_darcy*N[2];
const double crhs41 = N[2]*crhs20;
const double crhs42 = rho*(DN(2,0)*crhs8 + DN(2,1)*crhs9 + DN(2,2)*crhs10);
const double crhs43 = K_darcy*N[3];
const double crhs44 = N[3]*crhs20;
const double crhs45 = rho*(DN(3,0)*crhs8 + DN(3,1)*crhs9 + DN(3,2)*crhs10);
rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs16 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] + N[0]*crhs1 - N[0]*crhs11 - N[0]*crhs3 - N[0]*crhs5 - N[0]*crhs6 + crhs18*crhs19 - crhs18*crhs21 - crhs18*crhs22;
rhs[1]=-DN(0,0)*stress[3] + DN(0,1)*crhs0 - DN(0,1)*crhs16 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + N[0]*crhs23 - N[0]*crhs25 - N[0]*crhs26 - N[0]*crhs27 - N[0]*crhs28 + crhs19*crhs29 - crhs21*crhs29 - crhs22*crhs29;
rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] + DN(0,2)*crhs0 - DN(0,2)*crhs16 - DN(0,2)*stress[2] + N[0]*crhs30 - N[0]*crhs32 - N[0]*crhs33 - N[0]*crhs34 - N[0]*crhs35 + crhs19*crhs36 - crhs21*crhs36 - crhs22*crhs36;
rhs[3]=-DN(0,0)*crhs18 - DN(0,1)*crhs29 - DN(0,2)*crhs36 - N[0]*crhs14;
rhs[4]=DN(1,0)*crhs0 - DN(1,0)*crhs16 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] + N[1]*crhs1 - N[1]*crhs11 - N[1]*crhs3 - N[1]*crhs5 - N[1]*crhs6 + crhs18*crhs37 - crhs18*crhs38 - crhs18*crhs39;
rhs[5]=-DN(1,0)*stress[3] + DN(1,1)*crhs0 - DN(1,1)*crhs16 - DN(1,1)*stress[1] - DN(1,2)*stress[4] + N[1]*crhs23 - N[1]*crhs25 - N[1]*crhs26 - N[1]*crhs27 - N[1]*crhs28 + crhs29*crhs37 - crhs29*crhs38 - crhs29*crhs39;
rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] + DN(1,2)*crhs0 - DN(1,2)*crhs16 - DN(1,2)*stress[2] + N[1]*crhs30 - N[1]*crhs32 - N[1]*crhs33 - N[1]*crhs34 - N[1]*crhs35 + crhs36*crhs37 - crhs36*crhs38 - crhs36*crhs39;
rhs[7]=-DN(1,0)*crhs18 - DN(1,1)*crhs29 - DN(1,2)*crhs36 - N[1]*crhs14;
rhs[8]=DN(2,0)*crhs0 - DN(2,0)*crhs16 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] + N[2]*crhs1 - N[2]*crhs11 - N[2]*crhs3 - N[2]*crhs5 - N[2]*crhs6 + crhs18*crhs40 - crhs18*crhs41 - crhs18*crhs42;
rhs[9]=-DN(2,0)*stress[3] + DN(2,1)*crhs0 - DN(2,1)*crhs16 - DN(2,1)*stress[1] - DN(2,2)*stress[4] + N[2]*crhs23 - N[2]*crhs25 - N[2]*crhs26 - N[2]*crhs27 - N[2]*crhs28 + crhs29*crhs40 - crhs29*crhs41 - crhs29*crhs42;
rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] + DN(2,2)*crhs0 - DN(2,2)*crhs16 - DN(2,2)*stress[2] + N[2]*crhs30 - N[2]*crhs32 - N[2]*crhs33 - N[2]*crhs34 - N[2]*crhs35 + crhs36*crhs40 - crhs36*crhs41 - crhs36*crhs42;
rhs[11]=-DN(2,0)*crhs18 - DN(2,1)*crhs29 - DN(2,2)*crhs36 - N[2]*crhs14;
rhs[12]=DN(3,0)*crhs0 - DN(3,0)*crhs16 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] + N[3]*crhs1 - N[3]*crhs11 - N[3]*crhs3 - N[3]*crhs5 - N[3]*crhs6 + crhs18*crhs43 - crhs18*crhs44 - crhs18*crhs45;
rhs[13]=-DN(3,0)*stress[3] + DN(3,1)*crhs0 - DN(3,1)*crhs16 - DN(3,1)*stress[1] - DN(3,2)*stress[4] + N[3]*crhs23 - N[3]*crhs25 - N[3]*crhs26 - N[3]*crhs27 - N[3]*crhs28 + crhs29*crhs43 - crhs29*crhs44 - crhs29*crhs45;
rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] + DN(3,2)*crhs0 - DN(3,2)*crhs16 - DN(3,2)*stress[2] + N[3]*crhs30 - N[3]*crhs32 - N[3]*crhs33 - N[3]*crhs34 - N[3]*crhs35 + crhs36*crhs43 - crhs36*crhs44 - crhs36*crhs45;
rhs[15]=-DN(3,0)*crhs18 - DN(3,1)*crhs29 - DN(3,2)*crhs36 - N[3]*crhs14;


    noalias(rRHS) += rData.Weight * rhs;
}

template <>
void TwoFluidNavierStokes<TwoFluidNavierStokesData<2, 3>>::ComputeGaussPointEnrichmentContributions(
    TwoFluidNavierStokesData<2, 3> &rData,
    MatrixType &rV,
    MatrixType &rH,
    MatrixType &rKee,
    VectorType &rRHS_ee)
{

    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double bdf1 = rData.bdf1;
    const double bdf2 = rData.bdf2;

    const double dyn_tau = rData.DynamicTau;
    const double K_darcy = rData.DarcyTerm;

    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vnn = rData.Velocity_OldStep2;
    const auto &vmesh = rData.MeshVelocity;
    const auto &vconv = v - vmesh;
    const auto &f = rData.BodyForce;
    const auto &p = rData.Pressure;

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;
    const auto &Nenr = rData.Nenr;
    const auto &DNenr = rData.DN_DXenr;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    // Mass correction term
    const double volume_error_ratio = rData.VolumeError;

    auto &V = rData.V;
    auto &H = rData.H;
    auto &Kee = rData.Kee;
    auto &rhs_ee = rData.rhs_ee;

    array_1d<double, NumNodes> penr = ZeroVector(NumNodes); //penriched is considered to be zero as we do not want to store it

    const double cV0 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double cV1 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double cV2 = 1.0/(K_darcy + rho*stab_c2*sqrt(pow(cV0, 2) + pow(cV1, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cV3 = DNenr(0,0)*cV2;
const double cV4 = K_darcy*N[0];
const double cV5 = rho*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1));
const double cV6 = N[0]*cV5;
const double cV7 = rho*(DN(0,0)*cV0 + DN(0,1)*cV1);
const double cV8 = DNenr(1,0)*cV2;
const double cV9 = DNenr(2,0)*cV2;
const double cV10 = DNenr(0,1)*cV2;
const double cV11 = DNenr(1,1)*cV2;
const double cV12 = DNenr(2,1)*cV2;
const double cV13 = K_darcy*N[1];
const double cV14 = N[1]*cV5;
const double cV15 = rho*(DN(1,0)*cV0 + DN(1,1)*cV1);
const double cV16 = K_darcy*N[2];
const double cV17 = N[2]*cV5;
const double cV18 = rho*(DN(2,0)*cV0 + DN(2,1)*cV1);
V(0,0)=-DN(0,0)*Nenr[0] - cV3*cV4 + cV3*cV6 + cV3*cV7;
V(0,1)=-DN(0,0)*Nenr[1] - cV4*cV8 + cV6*cV8 + cV7*cV8;
V(0,2)=-DN(0,0)*Nenr[2] - cV4*cV9 + cV6*cV9 + cV7*cV9;
V(1,0)=-DN(0,1)*Nenr[0] - cV10*cV4 + cV10*cV6 + cV10*cV7;
V(1,1)=-DN(0,1)*Nenr[1] - cV11*cV4 + cV11*cV6 + cV11*cV7;
V(1,2)=-DN(0,1)*Nenr[2] - cV12*cV4 + cV12*cV6 + cV12*cV7;
V(2,0)=cV2*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1));
V(2,1)=cV2*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1));
V(2,2)=cV2*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1));
V(3,0)=-DN(1,0)*Nenr[0] - cV13*cV3 + cV14*cV3 + cV15*cV3;
V(3,1)=-DN(1,0)*Nenr[1] - cV13*cV8 + cV14*cV8 + cV15*cV8;
V(3,2)=-DN(1,0)*Nenr[2] - cV13*cV9 + cV14*cV9 + cV15*cV9;
V(4,0)=-DN(1,1)*Nenr[0] - cV10*cV13 + cV10*cV14 + cV10*cV15;
V(4,1)=-DN(1,1)*Nenr[1] - cV11*cV13 + cV11*cV14 + cV11*cV15;
V(4,2)=-DN(1,1)*Nenr[2] - cV12*cV13 + cV12*cV14 + cV12*cV15;
V(5,0)=cV2*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1));
V(5,1)=cV2*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1));
V(5,2)=cV2*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1));
V(6,0)=-DN(2,0)*Nenr[0] - cV16*cV3 + cV17*cV3 + cV18*cV3;
V(6,1)=-DN(2,0)*Nenr[1] - cV16*cV8 + cV17*cV8 + cV18*cV8;
V(6,2)=-DN(2,0)*Nenr[2] - cV16*cV9 + cV17*cV9 + cV18*cV9;
V(7,0)=-DN(2,1)*Nenr[0] - cV10*cV16 + cV10*cV17 + cV10*cV18;
V(7,1)=-DN(2,1)*Nenr[1] - cV11*cV16 + cV11*cV17 + cV11*cV18;
V(7,2)=-DN(2,1)*Nenr[2] - cV12*cV16 + cV12*cV17 + cV12*cV18;
V(8,0)=cV2*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1));
V(8,1)=cV2*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1));
V(8,2)=cV2*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1));


    const double cH0 = rho*volume_error_ratio;
const double cH1 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double cH2 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double cH3 = 1.0/(K_darcy + rho*stab_c2*sqrt(pow(cH1, 2) + pow(cH2, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cH4 = cH3*(K_darcy*N[0] + N[0]*cH0 + rho*(DN(0,0)*cH1 + DN(0,1)*cH2 + N[0]*bdf0));
const double cH5 = cH3*(K_darcy*N[1] + N[1]*cH0 + rho*(DN(1,0)*cH1 + DN(1,1)*cH2 + N[1]*bdf0));
const double cH6 = cH3*(K_darcy*N[2] + N[2]*cH0 + rho*(DN(2,0)*cH1 + DN(2,1)*cH2 + N[2]*bdf0));
H(0,0)=DN(0,0)*Nenr[0] + DNenr(0,0)*cH4;
H(0,1)=DN(0,1)*Nenr[0] + DNenr(0,1)*cH4;
H(0,2)=cH3*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1));
H(0,3)=DN(1,0)*Nenr[0] + DNenr(0,0)*cH5;
H(0,4)=DN(1,1)*Nenr[0] + DNenr(0,1)*cH5;
H(0,5)=cH3*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1));
H(0,6)=DN(2,0)*Nenr[0] + DNenr(0,0)*cH6;
H(0,7)=DN(2,1)*Nenr[0] + DNenr(0,1)*cH6;
H(0,8)=cH3*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1));
H(1,0)=DN(0,0)*Nenr[1] + DNenr(1,0)*cH4;
H(1,1)=DN(0,1)*Nenr[1] + DNenr(1,1)*cH4;
H(1,2)=cH3*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1));
H(1,3)=DN(1,0)*Nenr[1] + DNenr(1,0)*cH5;
H(1,4)=DN(1,1)*Nenr[1] + DNenr(1,1)*cH5;
H(1,5)=cH3*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1));
H(1,6)=DN(2,0)*Nenr[1] + DNenr(1,0)*cH6;
H(1,7)=DN(2,1)*Nenr[1] + DNenr(1,1)*cH6;
H(1,8)=cH3*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1));
H(2,0)=DN(0,0)*Nenr[2] + DNenr(2,0)*cH4;
H(2,1)=DN(0,1)*Nenr[2] + DNenr(2,1)*cH4;
H(2,2)=cH3*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1));
H(2,3)=DN(1,0)*Nenr[2] + DNenr(2,0)*cH5;
H(2,4)=DN(1,1)*Nenr[2] + DNenr(2,1)*cH5;
H(2,5)=cH3*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1));
H(2,6)=DN(2,0)*Nenr[2] + DNenr(2,0)*cH6;
H(2,7)=DN(2,1)*Nenr[2] + DNenr(2,1)*cH6;
H(2,8)=cH3*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1));


    const double cKee0 = 1.0/(K_darcy + rho*stab_c2*sqrt(pow(N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0), 2) + pow(N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1), 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cKee1 = cKee0*(DNenr(0,0)*DNenr(1,0) + DNenr(0,1)*DNenr(1,1));
const double cKee2 = cKee0*(DNenr(0,0)*DNenr(2,0) + DNenr(0,1)*DNenr(2,1));
const double cKee3 = cKee0*(DNenr(1,0)*DNenr(2,0) + DNenr(1,1)*DNenr(2,1));
Kee(0,0)=cKee0*(pow(DNenr(0,0), 2) + pow(DNenr(0,1), 2));
Kee(0,1)=cKee1;
Kee(0,2)=cKee2;
Kee(1,0)=cKee1;
Kee(1,1)=cKee0*(pow(DNenr(1,0), 2) + pow(DNenr(1,1), 2));
Kee(1,2)=cKee3;
Kee(2,0)=cKee2;
Kee(2,1)=cKee3;
Kee(2,2)=cKee0*(pow(DNenr(2,0), 2) + pow(DNenr(2,1), 2));


    const double crhs_ee0 = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
const double crhs_ee1 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
const double crhs_ee2 = crhs_ee0 + crhs_ee1 - volume_error_ratio;
const double crhs_ee3 = N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0);
const double crhs_ee4 = rho*volume_error_ratio;
const double crhs_ee5 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs_ee6 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs_ee7 = 1.0/(K_darcy + rho*stab_c2*sqrt(pow(crhs_ee5, 2) + pow(crhs_ee6, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs_ee8 = crhs_ee7*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DNenr(0,0)*penr[0] + DNenr(1,0)*penr[1] + DNenr(2,0)*penr[2] + K_darcy*crhs_ee3 + crhs_ee3*crhs_ee4 - rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0)) + rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + crhs_ee0*crhs_ee5 + crhs_ee6*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0))));
const double crhs_ee9 = N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1);
const double crhs_ee10 = crhs_ee7*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DNenr(0,1)*penr[0] + DNenr(1,1)*penr[1] + DNenr(2,1)*penr[2] + K_darcy*crhs_ee9 + crhs_ee4*crhs_ee9 - rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1)) + rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + crhs_ee1*crhs_ee6 + crhs_ee5*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1))));
rhs_ee[0]=-DNenr(0,0)*crhs_ee8 - DNenr(0,1)*crhs_ee10 - Nenr[0]*crhs_ee2;
rhs_ee[1]=-DNenr(1,0)*crhs_ee8 - DNenr(1,1)*crhs_ee10 - Nenr[1]*crhs_ee2;
rhs_ee[2]=-DNenr(2,0)*crhs_ee8 - DNenr(2,1)*crhs_ee10 - Nenr[2]*crhs_ee2;


    noalias(rV) += rData.Weight * V;
    noalias(rH) += rData.Weight * H;
    noalias(rKee) += rData.Weight * Kee;
    noalias(rRHS_ee) += rData.Weight * rhs_ee;
}

template <>
void TwoFluidNavierStokes<TwoFluidNavierStokesData<3, 4>>::ComputeGaussPointEnrichmentContributions(
    TwoFluidNavierStokesData<3, 4> &rData,
    MatrixType &rV,
    MatrixType &rH,
    MatrixType &rKee,
    VectorType &rRHS_ee)
{

    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double bdf1 = rData.bdf1;
    const double bdf2 = rData.bdf2;

    const double dyn_tau = rData.DynamicTau;
    const double K_darcy = rData.DarcyTerm;

    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vnn = rData.Velocity_OldStep2;
    const auto &vmesh = rData.MeshVelocity;
    const auto &vconv = v - vmesh;
    const auto &f = rData.BodyForce;
    const auto &p = rData.Pressure;
    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;
    const auto &Nenr = rData.Nenr;
    const auto &DNenr = rData.DN_DXenr;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    // Mass correction term
    const double volume_error_ratio = rData.VolumeError;

    auto &V = rData.V;
    auto &H = rData.H;
    auto &Kee = rData.Kee;
    auto &rhs_ee = rData.rhs_ee;

    array_1d<double, NumNodes> penr = ZeroVector(NumNodes); //penriched is considered to be zero as we do not want to store it

    const double cV0 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double cV1 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double cV2 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double cV3 = 1.0/(K_darcy + rho*stab_c2*sqrt(pow(cV0, 2) + pow(cV1, 2) + pow(cV2, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cV4 = DNenr(0,0)*cV3;
const double cV5 = K_darcy*N[0];
const double cV6 = rho*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2));
const double cV7 = N[0]*cV6;
const double cV8 = rho*(DN(0,0)*cV0 + DN(0,1)*cV1 + DN(0,2)*cV2);
const double cV9 = DNenr(1,0)*cV3;
const double cV10 = DNenr(2,0)*cV3;
const double cV11 = DNenr(3,0)*cV3;
const double cV12 = DNenr(0,1)*cV3;
const double cV13 = DNenr(1,1)*cV3;
const double cV14 = DNenr(2,1)*cV3;
const double cV15 = DNenr(3,1)*cV3;
const double cV16 = DNenr(0,2)*cV3;
const double cV17 = DNenr(1,2)*cV3;
const double cV18 = DNenr(2,2)*cV3;
const double cV19 = DNenr(3,2)*cV3;
const double cV20 = K_darcy*N[1];
const double cV21 = N[1]*cV6;
const double cV22 = rho*(DN(1,0)*cV0 + DN(1,1)*cV1 + DN(1,2)*cV2);
const double cV23 = K_darcy*N[2];
const double cV24 = N[2]*cV6;
const double cV25 = rho*(DN(2,0)*cV0 + DN(2,1)*cV1 + DN(2,2)*cV2);
const double cV26 = K_darcy*N[3];
const double cV27 = N[3]*cV6;
const double cV28 = rho*(DN(3,0)*cV0 + DN(3,1)*cV1 + DN(3,2)*cV2);
V(0,0)=-DN(0,0)*Nenr[0] - cV4*cV5 + cV4*cV7 + cV4*cV8;
V(0,1)=-DN(0,0)*Nenr[1] - cV5*cV9 + cV7*cV9 + cV8*cV9;
V(0,2)=-DN(0,0)*Nenr[2] - cV10*cV5 + cV10*cV7 + cV10*cV8;
V(0,3)=-DN(0,0)*Nenr[3] - cV11*cV5 + cV11*cV7 + cV11*cV8;
V(1,0)=-DN(0,1)*Nenr[0] - cV12*cV5 + cV12*cV7 + cV12*cV8;
V(1,1)=-DN(0,1)*Nenr[1] - cV13*cV5 + cV13*cV7 + cV13*cV8;
V(1,2)=-DN(0,1)*Nenr[2] - cV14*cV5 + cV14*cV7 + cV14*cV8;
V(1,3)=-DN(0,1)*Nenr[3] - cV15*cV5 + cV15*cV7 + cV15*cV8;
V(2,0)=-DN(0,2)*Nenr[0] - cV16*cV5 + cV16*cV7 + cV16*cV8;
V(2,1)=-DN(0,2)*Nenr[1] - cV17*cV5 + cV17*cV7 + cV17*cV8;
V(2,2)=-DN(0,2)*Nenr[2] - cV18*cV5 + cV18*cV7 + cV18*cV8;
V(2,3)=-DN(0,2)*Nenr[3] - cV19*cV5 + cV19*cV7 + cV19*cV8;
V(3,0)=cV3*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1) + DN(0,2)*DNenr(0,2));
V(3,1)=cV3*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1) + DN(0,2)*DNenr(1,2));
V(3,2)=cV3*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1) + DN(0,2)*DNenr(2,2));
V(3,3)=cV3*(DN(0,0)*DNenr(3,0) + DN(0,1)*DNenr(3,1) + DN(0,2)*DNenr(3,2));
V(4,0)=-DN(1,0)*Nenr[0] - cV20*cV4 + cV21*cV4 + cV22*cV4;
V(4,1)=-DN(1,0)*Nenr[1] - cV20*cV9 + cV21*cV9 + cV22*cV9;
V(4,2)=-DN(1,0)*Nenr[2] - cV10*cV20 + cV10*cV21 + cV10*cV22;
V(4,3)=-DN(1,0)*Nenr[3] - cV11*cV20 + cV11*cV21 + cV11*cV22;
V(5,0)=-DN(1,1)*Nenr[0] - cV12*cV20 + cV12*cV21 + cV12*cV22;
V(5,1)=-DN(1,1)*Nenr[1] - cV13*cV20 + cV13*cV21 + cV13*cV22;
V(5,2)=-DN(1,1)*Nenr[2] - cV14*cV20 + cV14*cV21 + cV14*cV22;
V(5,3)=-DN(1,1)*Nenr[3] - cV15*cV20 + cV15*cV21 + cV15*cV22;
V(6,0)=-DN(1,2)*Nenr[0] - cV16*cV20 + cV16*cV21 + cV16*cV22;
V(6,1)=-DN(1,2)*Nenr[1] - cV17*cV20 + cV17*cV21 + cV17*cV22;
V(6,2)=-DN(1,2)*Nenr[2] - cV18*cV20 + cV18*cV21 + cV18*cV22;
V(6,3)=-DN(1,2)*Nenr[3] - cV19*cV20 + cV19*cV21 + cV19*cV22;
V(7,0)=cV3*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1) + DN(1,2)*DNenr(0,2));
V(7,1)=cV3*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1) + DN(1,2)*DNenr(1,2));
V(7,2)=cV3*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1) + DN(1,2)*DNenr(2,2));
V(7,3)=cV3*(DN(1,0)*DNenr(3,0) + DN(1,1)*DNenr(3,1) + DN(1,2)*DNenr(3,2));
V(8,0)=-DN(2,0)*Nenr[0] - cV23*cV4 + cV24*cV4 + cV25*cV4;
V(8,1)=-DN(2,0)*Nenr[1] - cV23*cV9 + cV24*cV9 + cV25*cV9;
V(8,2)=-DN(2,0)*Nenr[2] - cV10*cV23 + cV10*cV24 + cV10*cV25;
V(8,3)=-DN(2,0)*Nenr[3] - cV11*cV23 + cV11*cV24 + cV11*cV25;
V(9,0)=-DN(2,1)*Nenr[0] - cV12*cV23 + cV12*cV24 + cV12*cV25;
V(9,1)=-DN(2,1)*Nenr[1] - cV13*cV23 + cV13*cV24 + cV13*cV25;
V(9,2)=-DN(2,1)*Nenr[2] - cV14*cV23 + cV14*cV24 + cV14*cV25;
V(9,3)=-DN(2,1)*Nenr[3] - cV15*cV23 + cV15*cV24 + cV15*cV25;
V(10,0)=-DN(2,2)*Nenr[0] - cV16*cV23 + cV16*cV24 + cV16*cV25;
V(10,1)=-DN(2,2)*Nenr[1] - cV17*cV23 + cV17*cV24 + cV17*cV25;
V(10,2)=-DN(2,2)*Nenr[2] - cV18*cV23 + cV18*cV24 + cV18*cV25;
V(10,3)=-DN(2,2)*Nenr[3] - cV19*cV23 + cV19*cV24 + cV19*cV25;
V(11,0)=cV3*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1) + DN(2,2)*DNenr(0,2));
V(11,1)=cV3*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1) + DN(2,2)*DNenr(1,2));
V(11,2)=cV3*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1) + DN(2,2)*DNenr(2,2));
V(11,3)=cV3*(DN(2,0)*DNenr(3,0) + DN(2,1)*DNenr(3,1) + DN(2,2)*DNenr(3,2));
V(12,0)=-DN(3,0)*Nenr[0] - cV26*cV4 + cV27*cV4 + cV28*cV4;
V(12,1)=-DN(3,0)*Nenr[1] - cV26*cV9 + cV27*cV9 + cV28*cV9;
V(12,2)=-DN(3,0)*Nenr[2] - cV10*cV26 + cV10*cV27 + cV10*cV28;
V(12,3)=-DN(3,0)*Nenr[3] - cV11*cV26 + cV11*cV27 + cV11*cV28;
V(13,0)=-DN(3,1)*Nenr[0] - cV12*cV26 + cV12*cV27 + cV12*cV28;
V(13,1)=-DN(3,1)*Nenr[1] - cV13*cV26 + cV13*cV27 + cV13*cV28;
V(13,2)=-DN(3,1)*Nenr[2] - cV14*cV26 + cV14*cV27 + cV14*cV28;
V(13,3)=-DN(3,1)*Nenr[3] - cV15*cV26 + cV15*cV27 + cV15*cV28;
V(14,0)=-DN(3,2)*Nenr[0] - cV16*cV26 + cV16*cV27 + cV16*cV28;
V(14,1)=-DN(3,2)*Nenr[1] - cV17*cV26 + cV17*cV27 + cV17*cV28;
V(14,2)=-DN(3,2)*Nenr[2] - cV18*cV26 + cV18*cV27 + cV18*cV28;
V(14,3)=-DN(3,2)*Nenr[3] - cV19*cV26 + cV19*cV27 + cV19*cV28;
V(15,0)=cV3*(DN(3,0)*DNenr(0,0) + DN(3,1)*DNenr(0,1) + DN(3,2)*DNenr(0,2));
V(15,1)=cV3*(DN(3,0)*DNenr(1,0) + DN(3,1)*DNenr(1,1) + DN(3,2)*DNenr(1,2));
V(15,2)=cV3*(DN(3,0)*DNenr(2,0) + DN(3,1)*DNenr(2,1) + DN(3,2)*DNenr(2,2));
V(15,3)=cV3*(DN(3,0)*DNenr(3,0) + DN(3,1)*DNenr(3,1) + DN(3,2)*DNenr(3,2));


    const double cH0 = rho*volume_error_ratio;
const double cH1 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double cH2 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double cH3 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double cH4 = 1.0/(K_darcy + rho*stab_c2*sqrt(pow(cH1, 2) + pow(cH2, 2) + pow(cH3, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cH5 = cH4*(K_darcy*N[0] + N[0]*cH0 + rho*(DN(0,0)*cH1 + DN(0,1)*cH2 + DN(0,2)*cH3 + N[0]*bdf0));
const double cH6 = cH4*(K_darcy*N[1] + N[1]*cH0 + rho*(DN(1,0)*cH1 + DN(1,1)*cH2 + DN(1,2)*cH3 + N[1]*bdf0));
const double cH7 = cH4*(K_darcy*N[2] + N[2]*cH0 + rho*(DN(2,0)*cH1 + DN(2,1)*cH2 + DN(2,2)*cH3 + N[2]*bdf0));
const double cH8 = cH4*(K_darcy*N[3] + N[3]*cH0 + rho*(DN(3,0)*cH1 + DN(3,1)*cH2 + DN(3,2)*cH3 + N[3]*bdf0));
H(0,0)=DN(0,0)*Nenr[0] + DNenr(0,0)*cH5;
H(0,1)=DN(0,1)*Nenr[0] + DNenr(0,1)*cH5;
H(0,2)=DN(0,2)*Nenr[0] + DNenr(0,2)*cH5;
H(0,3)=cH4*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1) + DN(0,2)*DNenr(0,2));
H(0,4)=DN(1,0)*Nenr[0] + DNenr(0,0)*cH6;
H(0,5)=DN(1,1)*Nenr[0] + DNenr(0,1)*cH6;
H(0,6)=DN(1,2)*Nenr[0] + DNenr(0,2)*cH6;
H(0,7)=cH4*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1) + DN(1,2)*DNenr(0,2));
H(0,8)=DN(2,0)*Nenr[0] + DNenr(0,0)*cH7;
H(0,9)=DN(2,1)*Nenr[0] + DNenr(0,1)*cH7;
H(0,10)=DN(2,2)*Nenr[0] + DNenr(0,2)*cH7;
H(0,11)=cH4*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1) + DN(2,2)*DNenr(0,2));
H(0,12)=DN(3,0)*Nenr[0] + DNenr(0,0)*cH8;
H(0,13)=DN(3,1)*Nenr[0] + DNenr(0,1)*cH8;
H(0,14)=DN(3,2)*Nenr[0] + DNenr(0,2)*cH8;
H(0,15)=cH4*(DN(3,0)*DNenr(0,0) + DN(3,1)*DNenr(0,1) + DN(3,2)*DNenr(0,2));
H(1,0)=DN(0,0)*Nenr[1] + DNenr(1,0)*cH5;
H(1,1)=DN(0,1)*Nenr[1] + DNenr(1,1)*cH5;
H(1,2)=DN(0,2)*Nenr[1] + DNenr(1,2)*cH5;
H(1,3)=cH4*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1) + DN(0,2)*DNenr(1,2));
H(1,4)=DN(1,0)*Nenr[1] + DNenr(1,0)*cH6;
H(1,5)=DN(1,1)*Nenr[1] + DNenr(1,1)*cH6;
H(1,6)=DN(1,2)*Nenr[1] + DNenr(1,2)*cH6;
H(1,7)=cH4*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1) + DN(1,2)*DNenr(1,2));
H(1,8)=DN(2,0)*Nenr[1] + DNenr(1,0)*cH7;
H(1,9)=DN(2,1)*Nenr[1] + DNenr(1,1)*cH7;
H(1,10)=DN(2,2)*Nenr[1] + DNenr(1,2)*cH7;
H(1,11)=cH4*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1) + DN(2,2)*DNenr(1,2));
H(1,12)=DN(3,0)*Nenr[1] + DNenr(1,0)*cH8;
H(1,13)=DN(3,1)*Nenr[1] + DNenr(1,1)*cH8;
H(1,14)=DN(3,2)*Nenr[1] + DNenr(1,2)*cH8;
H(1,15)=cH4*(DN(3,0)*DNenr(1,0) + DN(3,1)*DNenr(1,1) + DN(3,2)*DNenr(1,2));
H(2,0)=DN(0,0)*Nenr[2] + DNenr(2,0)*cH5;
H(2,1)=DN(0,1)*Nenr[2] + DNenr(2,1)*cH5;
H(2,2)=DN(0,2)*Nenr[2] + DNenr(2,2)*cH5;
H(2,3)=cH4*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1) + DN(0,2)*DNenr(2,2));
H(2,4)=DN(1,0)*Nenr[2] + DNenr(2,0)*cH6;
H(2,5)=DN(1,1)*Nenr[2] + DNenr(2,1)*cH6;
H(2,6)=DN(1,2)*Nenr[2] + DNenr(2,2)*cH6;
H(2,7)=cH4*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1) + DN(1,2)*DNenr(2,2));
H(2,8)=DN(2,0)*Nenr[2] + DNenr(2,0)*cH7;
H(2,9)=DN(2,1)*Nenr[2] + DNenr(2,1)*cH7;
H(2,10)=DN(2,2)*Nenr[2] + DNenr(2,2)*cH7;
H(2,11)=cH4*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1) + DN(2,2)*DNenr(2,2));
H(2,12)=DN(3,0)*Nenr[2] + DNenr(2,0)*cH8;
H(2,13)=DN(3,1)*Nenr[2] + DNenr(2,1)*cH8;
H(2,14)=DN(3,2)*Nenr[2] + DNenr(2,2)*cH8;
H(2,15)=cH4*(DN(3,0)*DNenr(2,0) + DN(3,1)*DNenr(2,1) + DN(3,2)*DNenr(2,2));
H(3,0)=DN(0,0)*Nenr[3] + DNenr(3,0)*cH5;
H(3,1)=DN(0,1)*Nenr[3] + DNenr(3,1)*cH5;
H(3,2)=DN(0,2)*Nenr[3] + DNenr(3,2)*cH5;
H(3,3)=cH4*(DN(0,0)*DNenr(3,0) + DN(0,1)*DNenr(3,1) + DN(0,2)*DNenr(3,2));
H(3,4)=DN(1,0)*Nenr[3] + DNenr(3,0)*cH6;
H(3,5)=DN(1,1)*Nenr[3] + DNenr(3,1)*cH6;
H(3,6)=DN(1,2)*Nenr[3] + DNenr(3,2)*cH6;
H(3,7)=cH4*(DN(1,0)*DNenr(3,0) + DN(1,1)*DNenr(3,1) + DN(1,2)*DNenr(3,2));
H(3,8)=DN(2,0)*Nenr[3] + DNenr(3,0)*cH7;
H(3,9)=DN(2,1)*Nenr[3] + DNenr(3,1)*cH7;
H(3,10)=DN(2,2)*Nenr[3] + DNenr(3,2)*cH7;
H(3,11)=cH4*(DN(2,0)*DNenr(3,0) + DN(2,1)*DNenr(3,1) + DN(2,2)*DNenr(3,2));
H(3,12)=DN(3,0)*Nenr[3] + DNenr(3,0)*cH8;
H(3,13)=DN(3,1)*Nenr[3] + DNenr(3,1)*cH8;
H(3,14)=DN(3,2)*Nenr[3] + DNenr(3,2)*cH8;
H(3,15)=cH4*(DN(3,0)*DNenr(3,0) + DN(3,1)*DNenr(3,1) + DN(3,2)*DNenr(3,2));


    const double cKee0 = 1.0/(K_darcy + rho*stab_c2*sqrt(pow(N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0), 2) + pow(N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1), 2) + pow(N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2), 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cKee1 = cKee0*(DNenr(0,0)*DNenr(1,0) + DNenr(0,1)*DNenr(1,1) + DNenr(0,2)*DNenr(1,2));
const double cKee2 = cKee0*(DNenr(0,0)*DNenr(2,0) + DNenr(0,1)*DNenr(2,1) + DNenr(0,2)*DNenr(2,2));
const double cKee3 = cKee0*(DNenr(0,0)*DNenr(3,0) + DNenr(0,1)*DNenr(3,1) + DNenr(0,2)*DNenr(3,2));
const double cKee4 = cKee0*(DNenr(1,0)*DNenr(2,0) + DNenr(1,1)*DNenr(2,1) + DNenr(1,2)*DNenr(2,2));
const double cKee5 = cKee0*(DNenr(1,0)*DNenr(3,0) + DNenr(1,1)*DNenr(3,1) + DNenr(1,2)*DNenr(3,2));
const double cKee6 = cKee0*(DNenr(2,0)*DNenr(3,0) + DNenr(2,1)*DNenr(3,1) + DNenr(2,2)*DNenr(3,2));
Kee(0,0)=cKee0*(pow(DNenr(0,0), 2) + pow(DNenr(0,1), 2) + pow(DNenr(0,2), 2));
Kee(0,1)=cKee1;
Kee(0,2)=cKee2;
Kee(0,3)=cKee3;
Kee(1,0)=cKee1;
Kee(1,1)=cKee0*(pow(DNenr(1,0), 2) + pow(DNenr(1,1), 2) + pow(DNenr(1,2), 2));
Kee(1,2)=cKee4;
Kee(1,3)=cKee5;
Kee(2,0)=cKee2;
Kee(2,1)=cKee4;
Kee(2,2)=cKee0*(pow(DNenr(2,0), 2) + pow(DNenr(2,1), 2) + pow(DNenr(2,2), 2));
Kee(2,3)=cKee6;
Kee(3,0)=cKee3;
Kee(3,1)=cKee5;
Kee(3,2)=cKee6;
Kee(3,3)=cKee0*(pow(DNenr(3,0), 2) + pow(DNenr(3,1), 2) + pow(DNenr(3,2), 2));


    const double crhs_ee0 = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0);
const double crhs_ee1 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
const double crhs_ee2 = DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
const double crhs_ee3 = crhs_ee0 + crhs_ee1 + crhs_ee2 - volume_error_ratio;
const double crhs_ee4 = N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0);
const double crhs_ee5 = rho*volume_error_ratio;
const double crhs_ee6 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs_ee7 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs_ee8 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs_ee9 = 1.0/(K_darcy + rho*stab_c2*sqrt(pow(crhs_ee6, 2) + pow(crhs_ee7, 2) + pow(crhs_ee8, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs_ee10 = crhs_ee9*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] + DNenr(0,0)*penr[0] + DNenr(1,0)*penr[1] + DNenr(2,0)*penr[2] + DNenr(3,0)*penr[3] + K_darcy*crhs_ee4 + crhs_ee4*crhs_ee5 - rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0)) + rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0)) + crhs_ee0*crhs_ee6 + crhs_ee7*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)) + crhs_ee8*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0))));
const double crhs_ee11 = N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1);
const double crhs_ee12 = crhs_ee9*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] + DNenr(0,1)*penr[0] + DNenr(1,1)*penr[1] + DNenr(2,1)*penr[2] + DNenr(3,1)*penr[3] + K_darcy*crhs_ee11 + crhs_ee11*crhs_ee5 - rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1)) + rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1)) + crhs_ee1*crhs_ee7 + crhs_ee6*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + crhs_ee8*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1))));
const double crhs_ee13 = N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2);
const double crhs_ee14 = crhs_ee9*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] + DNenr(0,2)*penr[0] + DNenr(1,2)*penr[1] + DNenr(2,2)*penr[2] + DNenr(3,2)*penr[3] + K_darcy*crhs_ee13 + crhs_ee13*crhs_ee5 - rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2)) + rho*(N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2)) + crhs_ee2*crhs_ee8 + crhs_ee6*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + crhs_ee7*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2))));
rhs_ee[0]=-DNenr(0,0)*crhs_ee10 - DNenr(0,1)*crhs_ee12 - DNenr(0,2)*crhs_ee14 - Nenr[0]*crhs_ee3;
rhs_ee[1]=-DNenr(1,0)*crhs_ee10 - DNenr(1,1)*crhs_ee12 - DNenr(1,2)*crhs_ee14 - Nenr[1]*crhs_ee3;
rhs_ee[2]=-DNenr(2,0)*crhs_ee10 - DNenr(2,1)*crhs_ee12 - DNenr(2,2)*crhs_ee14 - Nenr[2]*crhs_ee3;
rhs_ee[3]=-DNenr(3,0)*crhs_ee10 - DNenr(3,1)*crhs_ee12 - DNenr(3,2)*crhs_ee14 - Nenr[3]*crhs_ee3;


    noalias(rV) += rData.Weight * V;
    noalias(rH) += rData.Weight * H;
    noalias(rKee) += rData.Weight * Kee;
    noalias(rRHS_ee) += rData.Weight * rhs_ee;
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::ComputeSplitting(
    TElementData &rData,
    MatrixType &rShapeFunctionsPos,
    MatrixType &rShapeFunctionsNeg,
    MatrixType &rEnrichedShapeFunctionsPos,
    MatrixType &rEnrichedShapeFunctionsNeg,
    GeometryType::ShapeFunctionsGradientsType &rShapeDerivativesPos,
    GeometryType::ShapeFunctionsGradientsType &rShapeDerivativesNeg,
    GeometryType::ShapeFunctionsGradientsType &rEnrichedShapeDerivativesPos,
    GeometryType::ShapeFunctionsGradientsType &rEnrichedShapeDerivativesNeg,
    ModifiedShapeFunctions::Pointer pModifiedShapeFunctions)
{
    // Set the positive and negative enrichment interpolation matrices
    // Note that the enrichment is constructed using the standard shape functions such that:
    // In the negative distance region, the enrichment functions correspondig to the negative
    // distance nodes are null and the positive distance nodes are equal to the standard shape
    // functions. On the contrary, for the positive distance region, the enrichment functions
    // corresponding to the positive distance nodes are null meanwhile the negative distance
    // nodes are equal to the standard. This yields a discontinuous enrichment space.
    Matrix enr_neg_interp = ZeroMatrix(NumNodes, NumNodes);
    Matrix enr_pos_interp = ZeroMatrix(NumNodes, NumNodes);

    for (unsigned int i = 0; i < NumNodes; ++i){
        if (rData.Distance[i] > 0.0){
            enr_neg_interp(i, i) = 1.0;
        } else{
            enr_pos_interp(i, i) = 1.0;
        }
    }

    // Call the positive side modified shape functions calculator
    pModifiedShapeFunctions->ComputePositiveSideShapeFunctionsAndGradientsValues(
        rShapeFunctionsPos,
        rShapeDerivativesPos,
        rData.w_gauss_pos_side,
        GeometryData::IntegrationMethod::GI_GAUSS_2);

    // Call the negative side modified shape functions calculator
    pModifiedShapeFunctions->ComputeNegativeSideShapeFunctionsAndGradientsValues(
        rShapeFunctionsNeg,
        rShapeDerivativesNeg,
        rData.w_gauss_neg_side,
        GeometryData::IntegrationMethod::GI_GAUSS_2);

    // Compute the enrichment shape function values using the enrichment interpolation matrices
    rEnrichedShapeFunctionsPos = prod(rShapeFunctionsPos, enr_pos_interp);
    rEnrichedShapeFunctionsNeg = prod(rShapeFunctionsNeg, enr_neg_interp);

    // Compute the enrichment shape function gradient values using the enrichment interpolation matrices
    rEnrichedShapeDerivativesPos = rShapeDerivativesPos;
    rEnrichedShapeDerivativesNeg = rShapeDerivativesNeg;

    for (unsigned int i = 0; i < rShapeDerivativesPos.size(); ++i){
        rEnrichedShapeDerivativesPos[i] = prod(enr_pos_interp, rShapeDerivativesPos[i]);
    }

    for (unsigned int i = 0; i < rShapeDerivativesNeg.size(); ++i){
        rEnrichedShapeDerivativesNeg[i] = prod(enr_neg_interp, rShapeDerivativesNeg[i]);
    }

    rData.NumberOfDivisions = (pModifiedShapeFunctions->pGetSplittingUtil())->mDivisionsNumber;
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::ComputeSplitInterface(
    const TElementData &rData,
    MatrixType& rInterfaceShapeFunctionNeg,
    MatrixType& rEnrInterfaceShapeFunctionPos,
    MatrixType& rEnrInterfaceShapeFunctionNeg,
    GeometryType::ShapeFunctionsGradientsType& rInterfaceShapeDerivativesNeg,
    Vector& rInterfaceWeightsNeg,
    std::vector<array_1d<double,3>>& rInterfaceNormalsNeg,
    ModifiedShapeFunctions::Pointer pModifiedShapeFunctions)
{
    Matrix enr_neg_interp = ZeroMatrix(NumNodes, NumNodes);
    Matrix enr_pos_interp = ZeroMatrix(NumNodes, NumNodes);

    for (unsigned int i = 0; i < NumNodes; ++i){
        if (rData.Distance[i] > 0.0){
            enr_neg_interp(i, i) = 1.0;
        } else{
            enr_pos_interp(i, i) = 1.0;
        }
    }

    // Call the Interface negative side shape functions calculator
    pModifiedShapeFunctions->ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
        rInterfaceShapeFunctionNeg,
        rInterfaceShapeDerivativesNeg,
        rInterfaceWeightsNeg,
        GeometryData::IntegrationMethod::GI_GAUSS_2);

    // Call the Interface negative side normal functions calculator
    pModifiedShapeFunctions->ComputeNegativeSideInterfaceAreaNormals(
        rInterfaceNormalsNeg,
        GeometryData::IntegrationMethod::GI_GAUSS_2);

    for (unsigned int gp = 0; gp < rInterfaceNormalsNeg.size(); ++gp){
        const double normal_norm = norm_2(rInterfaceNormalsNeg[gp]);
        rInterfaceNormalsNeg[gp] /= normal_norm;
    }

    // Compute the enrichment shape function values at the interface gauss points using the enrichment interpolation matrices
    rEnrInterfaceShapeFunctionPos = prod(rInterfaceShapeFunctionNeg, enr_pos_interp);
    rEnrInterfaceShapeFunctionNeg = prod(rInterfaceShapeFunctionNeg, enr_neg_interp);
}

template <>
ModifiedShapeFunctions::UniquePointer TwoFluidNavierStokes< TwoFluidNavierStokesData<2, 3> >::pGetModifiedShapeFunctionsUtility(
    const GeometryType::Pointer pGeometry,
    const Vector& rDistances)
{
    return Kratos::make_unique<Triangle2D3ModifiedShapeFunctions>(pGeometry, rDistances);
}

template <>
ModifiedShapeFunctions::UniquePointer TwoFluidNavierStokes< TwoFluidNavierStokesData<3, 4> >::pGetModifiedShapeFunctionsUtility(
        const GeometryType::Pointer pGeometry,
        const Vector& rDistances)
{
    return Kratos::make_unique<Tetrahedra3D4ModifiedShapeFunctions>(pGeometry, rDistances);
}

template <>
ModifiedShapeFunctions::UniquePointer TwoFluidNavierStokes< TwoFluidNavierStokesAlphaMethodData<2, 3> >::pGetModifiedShapeFunctionsUtility(
    const GeometryType::Pointer pGeometry,
    const Vector& rDistances)
{
    return Kratos::make_unique<Triangle2D3ModifiedShapeFunctions>(pGeometry, rDistances);
}

template <>
ModifiedShapeFunctions::UniquePointer TwoFluidNavierStokes< TwoFluidNavierStokesAlphaMethodData<3, 4> >::pGetModifiedShapeFunctionsUtility(
        const GeometryType::Pointer pGeometry,
        const Vector& rDistances)
{
    return Kratos::make_unique<Tetrahedra3D4ModifiedShapeFunctions>(pGeometry, rDistances);
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::CalculateCurvatureOnInterfaceGaussPoints(
        const Matrix& rInterfaceShapeFunctions,
        Vector& rInterfaceCurvature)
{
    const auto& r_geom = this->GetGeometry();
    const unsigned int n_gpt = rInterfaceShapeFunctions.size1();

    rInterfaceCurvature.resize(n_gpt, false);

    for (unsigned int gpt = 0; gpt < n_gpt; ++gpt){
        double curvature = 0.0;
        for (unsigned int i = 0; i < NumNodes; ++i){
            curvature += rInterfaceShapeFunctions(gpt,i) * r_geom[i].GetValue(CURVATURE);
        }
        rInterfaceCurvature[gpt] = curvature;
    }
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::SurfaceTension(
    const double SurfaceTensionCoefficient,
    const Vector& rCurvature,
    const Vector& rInterfaceWeights,
    const Matrix& rInterfaceShapeFunctions,
    const std::vector<array_1d<double,3>>& rInterfaceNormalsNeg,
    VectorType& rRHS)
{
    for (unsigned int intgp = 0; intgp < rInterfaceWeights.size(); ++intgp){
        const double intgp_curv = rCurvature(intgp);
        const double intgp_w = rInterfaceWeights(intgp);
        const auto& intgp_normal = rInterfaceNormalsNeg[intgp];
        for (unsigned int j = 0; j < NumNodes; ++j){
            for (unsigned int dim = 0; dim < NumNodes-1; ++dim){
                rRHS[ j*(NumNodes) + dim ] -= SurfaceTensionCoefficient*intgp_normal[dim]
                    *intgp_curv*intgp_w*rInterfaceShapeFunctions(intgp,j);
            }
        }
    }
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::PressureGradientStabilization(
    const TElementData& rData,
    const Vector& rInterfaceWeights,
    const Matrix& rEnrInterfaceShapeFunctionPos,
    const Matrix& rEnrInterfaceShapeFunctionNeg,
    const GeometryType::ShapeFunctionsGradientsType& rInterfaceShapeDerivatives,
    MatrixType& rKeeTot,
    VectorType& rRHSeeTot)
{
    MatrixType kee = ZeroMatrix(NumNodes, NumNodes);
    VectorType rhs_enr = ZeroVector(NumNodes);

    Matrix enr_neg_interp = ZeroMatrix(NumNodes, NumNodes);
    Matrix enr_pos_interp = ZeroMatrix(NumNodes, NumNodes);

    double positive_density = 0.0;
    double negative_density = 0.0;
    double positive_viscosity = 0.0;
    double negative_viscosity = 0.0;

    for (unsigned int i = 0; i < NumNodes; ++i){
        if (rData.Distance[i] > 0.0){
            enr_neg_interp(i, i) = 1.0;
            positive_density = rData.NodalDensity[i];
            positive_viscosity = rData.NodalDynamicViscosity[i];
        } else{
            enr_pos_interp(i, i) = 1.0;
            negative_density = rData.NodalDensity[i];
            negative_viscosity = rData.NodalDynamicViscosity[i];
        }
    }

    GeometryType::ShapeFunctionsGradientsType EnrichedInterfaceShapeDerivativesPos = rInterfaceShapeDerivatives;
    GeometryType::ShapeFunctionsGradientsType EnrichedInterfaceShapeDerivativesNeg = rInterfaceShapeDerivatives;

    for (unsigned int i = 0; i < rInterfaceShapeDerivatives.size(); ++i){
        EnrichedInterfaceShapeDerivativesPos[i] = prod(enr_pos_interp, rInterfaceShapeDerivatives[i]);
    }

    for (unsigned int i = 0; i < rInterfaceShapeDerivatives.size(); ++i){
        EnrichedInterfaceShapeDerivativesNeg[i] = prod(enr_neg_interp, rInterfaceShapeDerivatives[i]);
    }

    double positive_volume = 0.0;
    double negative_volume = 0.0;
    for (unsigned int igauss_pos = 0; igauss_pos < rData.w_gauss_pos_side.size(); ++igauss_pos){
        positive_volume += rData.w_gauss_pos_side[igauss_pos];
    }

    for (unsigned int igauss_neg = 0; igauss_neg < rData.w_gauss_neg_side.size(); ++igauss_neg){
        negative_volume += rData.w_gauss_neg_side[igauss_neg];
    }
    const double element_volume = positive_volume + negative_volume;

    const auto& r_geom = this->GetGeometry();
    const double h_elem = rData.ElementSize;

    double cut_area = 0.0;
    for (unsigned int gp = 0; gp < rInterfaceWeights.size(); ++gp){
        cut_area += rInterfaceWeights[gp];
    }

    const double density = 1.0/(1.0/positive_density + 1.0/negative_density);
    const double viscosity = 1.0/(1.0/positive_viscosity + 1.0/negative_viscosity);

    // Stabilization parameters
    const double cut_stabilization_coefficient = 1.0;
    const double stab_c1 = 4.0;
    const double stab_c2 = 2.0;
    const double dyn_tau = rData.DynamicTau;

    const double dt = rData.DeltaTime;

    const auto v_convection = rData.Velocity - rData.MeshVelocity;

    for (unsigned int gp = 0; gp < rInterfaceWeights.size(); ++gp){

        Vector vconv = ZeroVector(Dim);
        double positive_weight = 0.0;
        double negative_weight = 0.0;

        for (unsigned int j = 0; j < NumNodes; ++j){
            for (unsigned int dim = 0; dim < Dim; ++dim){
                vconv[dim] += (rEnrInterfaceShapeFunctionNeg(gp, j) + rEnrInterfaceShapeFunctionPos(gp, j))
                    *v_convection(j,dim);
            }
            positive_weight += rEnrInterfaceShapeFunctionNeg(gp, j);
            negative_weight += rEnrInterfaceShapeFunctionPos(gp, j);
        }

        const double v_conv_norm = norm_2(vconv);

        const double penalty_coefficient = cut_stabilization_coefficient *
            density * 1.0 / (dyn_tau * density / dt + stab_c1 * viscosity / h_elem / h_elem +
                                stab_c2 * density * v_conv_norm / h_elem) * element_volume / cut_area;

        const auto& r_gp_enriched_interface_shape_derivatives_pos = EnrichedInterfaceShapeDerivativesPos[gp];
        const auto& r_gp_enriched_interface_shape_derivatives_neg = EnrichedInterfaceShapeDerivativesNeg[gp];

        for (unsigned int i = 0; i < NumNodes; ++i){

            for (unsigned int j = 0; j < NumNodes; ++j){

                const auto& r_pressure_gradient_j = r_geom[j].GetValue(PRESSURE_GRADIENT);

                for (unsigned int dim = 0; dim < Dim; ++dim){
                    kee(i, j) += penalty_coefficient * rInterfaceWeights[gp] *
                        ( r_gp_enriched_interface_shape_derivatives_pos(i,dim) - r_gp_enriched_interface_shape_derivatives_neg(i,dim) )*
                        ( r_gp_enriched_interface_shape_derivatives_pos(j,dim) - r_gp_enriched_interface_shape_derivatives_neg(j,dim) );

                    rhs_enr(i) += penalty_coefficient * rInterfaceWeights[gp] *
                        ( r_gp_enriched_interface_shape_derivatives_pos(i,dim) - r_gp_enriched_interface_shape_derivatives_neg(i,dim) )*
                        (rEnrInterfaceShapeFunctionNeg(gp, j)/positive_weight - rEnrInterfaceShapeFunctionPos(gp, j)/negative_weight)*
                        r_pressure_gradient_j(dim);
                }
            }
        }
    }

    noalias(rKeeTot) += kee;
    noalias(rRHSeeTot) += rhs_enr;
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::CalculateStrainRate(TElementData& rData) const
{
    FluidElement<TElementData>::CalculateStrainRate(rData);
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::CondenseEnrichmentWithContinuity(
    const TElementData &rData,
    Matrix &rLeftHandSideMatrix,
    VectorType &rRightHandSideVector,
    const MatrixType &rHtot,
    const MatrixType &rVtot,
    MatrixType &rKeeTot,
    const VectorType &rRHSeeTot)
{
    const double min_area_ratio = 1e-7;

    // Compute positive side, negative side and total volumes
    double positive_volume = 0.0;
    double negative_volume = 0.0;
    for (unsigned int igauss_pos = 0; igauss_pos < rData.w_gauss_pos_side.size(); ++igauss_pos){
        positive_volume += rData.w_gauss_pos_side[igauss_pos];
    }

    for (unsigned int igauss_neg = 0; igauss_neg < rData.w_gauss_neg_side.size(); ++igauss_neg){
        negative_volume += rData.w_gauss_neg_side[igauss_neg];
    }
    const double Vol = positive_volume + negative_volume;

    //We only enrich elements which are not almost empty/full
    if (positive_volume / Vol > min_area_ratio && negative_volume / Vol > min_area_ratio) {

        // Compute the maximum diagonal value in the enrichment stiffness matrix
        double max_diag = 0.0;
        for (unsigned int k = 0; k < NumNodes; ++k){
            if (std::abs(rKeeTot(k, k)) > max_diag){
                max_diag = std::abs(rKeeTot(k, k));
            }
        }
        if (max_diag == 0.0){
            max_diag = 1.0;
        }
        // "weakly" impose continuity
        for (unsigned int i = 0; i < Dim; ++i){
            const double di = std::abs(rData.Distance[i]);
            for (unsigned int j = i + 1; j < NumNodes; ++j){
                const double dj = std::abs(rData.Distance[j]);
                // Check if the edge is cut, if it is, set the penalty constraint
                if (rData.Distance[i] * rData.Distance[j] < 0.0){
                    double sum_d = di + dj;
                    double Ni = dj / sum_d;
                    double Nj = di / sum_d;
                    double penalty_coeff = max_diag * 0.001; // h/BDFVector[0];
                    rKeeTot(i, i) += penalty_coeff * Ni * Ni;
                    rKeeTot(i, j) -= penalty_coeff * Ni * Nj;
                    rKeeTot(j, i) -= penalty_coeff * Nj * Ni;
                    rKeeTot(j, j) += penalty_coeff * Nj * Nj;
                }
            }
        }

        // Enrichment condensation (add to LHS and RHS the enrichment contributions)
        double det;
        MatrixType inverse_diag(NumNodes, NumNodes);
        MathUtils<double>::InvertMatrix(rKeeTot, inverse_diag, det);

        const Matrix tmp = prod(inverse_diag, rHtot);
        noalias(rLeftHandSideMatrix) -= prod(rVtot, tmp);

        const Vector tmp2 = prod(inverse_diag, rRHSeeTot);
        noalias(rRightHandSideVector) -= prod(rVtot, tmp2);
    }
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::CondenseEnrichment(
    Matrix &rLeftHandSideMatrix,
    VectorType &rRightHandSideVector,
    const MatrixType &rHtot,
    const MatrixType &rVtot,
    MatrixType &rKeeTot,
    const VectorType &rRHSeeTot)
{
    // Enrichment condensation (add to LHS and RHS the enrichment contributions)
    double det;
    MatrixType inverse_diag(NumNodes, NumNodes);
    MathUtils<double>::InvertMatrix(rKeeTot, inverse_diag, det);

    const Matrix tmp = prod(inverse_diag, rHtot);
    noalias(rLeftHandSideMatrix) -= prod(rVtot, tmp);

    const Vector tmp2 = prod(inverse_diag, rRHSeeTot);
    noalias(rRightHandSideVector) -= prod(rVtot, tmp2);
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::AddSurfaceTensionContribution(
    const TElementData& rData,
    MatrixType& rInterfaceShapeFunction,
    MatrixType& rEnrInterfaceShapeFunctionPos,
    MatrixType& rEnrInterfaceShapeFunctionNeg,
    GeometryType::ShapeFunctionsGradientsType& rInterfaceShapeDerivatives,
    Vector& rInterfaceWeights,
    std::vector< array_1d<double, 3> >& rInterfaceNormalsNeg,
    Matrix &rLeftHandSideMatrix,
    VectorType &rRightHandSideVector,
    const MatrixType &rHtot,
    const MatrixType &rVtot,
    MatrixType &rKeeTot,
    VectorType &rRHSeeTot)
{
    // Surface tension coefficient is set in material properties
    const double surface_tension_coefficient = this->GetProperties().GetValue(SURFACE_TENSION_COEFFICIENT);

    Vector gauss_pts_curvature; // curvatures calculated on interface Gauss points

    CalculateCurvatureOnInterfaceGaussPoints(
        rInterfaceShapeFunction,
        gauss_pts_curvature);

    SurfaceTension(
        surface_tension_coefficient,
        gauss_pts_curvature,
        rInterfaceWeights,
        rInterfaceShapeFunction,
        rInterfaceNormalsNeg,
        rRightHandSideVector);

    this->PressureGradientStabilization(
        rData,
        rInterfaceWeights,
        rEnrInterfaceShapeFunctionPos,
        rEnrInterfaceShapeFunctionNeg,
        rInterfaceShapeDerivatives,
        rKeeTot,
        rRHSeeTot);

    CondenseEnrichment(rLeftHandSideMatrix, rRightHandSideVector, rHtot, rVtot, rKeeTot, rRHSeeTot);
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::save(Serializer &rSerializer) const
{
    using BaseType = FluidElement<TElementData>;
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::load(Serializer &rSerializer)
{
    using BaseType = FluidElement<TElementData>;
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
}


template <class TElementData>
void TwoFluidNavierStokes<TElementData>::CalculateOnIntegrationPoints(
    const Variable<double> &rVariable,
    std::vector<double> &rValues,
    const ProcessInfo &rCurrentProcessInfo )
{
    if (rVariable == DIVERGENCE){

        const auto& rGeom = this->GetGeometry();
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::IntegrationMethod::GI_GAUSS_2);
        const unsigned int num_gauss = IntegrationPoints.size();

        if (rValues.size() != num_gauss){
            rValues.resize(num_gauss);
        }

        Vector gauss_pts_jacobian_determinant = ZeroVector(num_gauss);
        GeometryData::ShapeFunctionsGradientsType DN_DX;
        rGeom.ShapeFunctionsIntegrationPointsGradients(DN_DX, gauss_pts_jacobian_determinant, GeometryData::IntegrationMethod::GI_GAUSS_2);

        for (unsigned int i_gauss = 0; i_gauss < num_gauss; ++i_gauss){

            const Matrix gp_DN_DX = DN_DX[i_gauss];
            double DVi_DXi = 0.0;

            for(unsigned int nnode = 0; nnode < NumNodes; ++nnode){

                const array_1d<double,3> vel = rGeom[nnode].GetSolutionStepValue(VELOCITY);
                for(unsigned int ndim = 0; ndim < Dim; ++ndim){
                    DVi_DXi += gp_DN_DX(nnode, ndim) * vel[ndim];
                }
            }
            rValues[i_gauss] = DVi_DXi;
        }
    }
}

template <>
void TwoFluidNavierStokes<TwoFluidNavierStokesAlphaMethodData<2, 3>>::ComputeGaussPointLHSContribution(
    TwoFluidNavierStokesAlphaMethodData<2, 3> &rData,
    MatrixType &rLHS)
{
    KRATOS_ERROR << "This function should never be called. It is only to enable the explicit template instantiation." << std::endl;
}

template <>
void TwoFluidNavierStokes<TwoFluidNavierStokesAlphaMethodData<3, 4>>::ComputeGaussPointLHSContribution(
    TwoFluidNavierStokesAlphaMethodData<3, 4> &rData,
    MatrixType &rLHS)
{
    KRATOS_ERROR << "This function should never be called. It is only to enable the explicit template instantiation." << std::endl;
}

template <>
void TwoFluidNavierStokes<TwoFluidNavierStokesAlphaMethodData<2, 3>>::ComputeGaussPointRHSContribution(
    TwoFluidNavierStokesAlphaMethodData<2, 3> &rData,
    VectorType &rRHS)
{
    KRATOS_ERROR << "This function should never be called. It is only to enable the explicit template instantiation." << std::endl;
}

template <>
void TwoFluidNavierStokes<TwoFluidNavierStokesAlphaMethodData<3, 4>>::ComputeGaussPointRHSContribution(
    TwoFluidNavierStokesAlphaMethodData<3, 4> &rData,
    VectorType &RLHS)
{
    KRATOS_ERROR << "This function should never be called. It is only to enable the explicit template instantiation." << std::endl;
}

template <>
void TwoFluidNavierStokes<TwoFluidNavierStokesAlphaMethodData<2, 3>>::ComputeGaussPointEnrichmentContributions(
    TwoFluidNavierStokesAlphaMethodData<2, 3> &rData,
    MatrixType &rV,
    MatrixType &rH,
    MatrixType &rKee,
    VectorType &rRHS_ee)
{
    KRATOS_ERROR << "This function should never be called. It is only to enable the explicit template instantiation." << std::endl;
}

template <>
void TwoFluidNavierStokes<TwoFluidNavierStokesAlphaMethodData<3, 4>>::ComputeGaussPointEnrichmentContributions(
    TwoFluidNavierStokesAlphaMethodData<3, 4> &rData,
    MatrixType &rV,
    MatrixType &rH,
    MatrixType &rKee,
    VectorType &rRHS_ee)
{
    KRATOS_ERROR << "This function should never be called. It is only to enable the explicit template instantiation." << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class TwoFluidNavierStokes<TwoFluidNavierStokesData<2, 3>>;
template class TwoFluidNavierStokes<TwoFluidNavierStokesData<3, 4>>;

template class TwoFluidNavierStokes<TwoFluidNavierStokesAlphaMethodData<2, 3>>;
template class TwoFluidNavierStokes<TwoFluidNavierStokesAlphaMethodData<3, 4>>;

} // namespace Kratos
