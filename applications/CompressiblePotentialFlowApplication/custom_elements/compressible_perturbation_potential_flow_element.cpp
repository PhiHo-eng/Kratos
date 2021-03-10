//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Inigo Lopez and Riccardo Rossi
//

#include "compressible_perturbation_potential_flow_element.h"
#include "compressible_potential_flow_application_variables.h"
#include "fluid_dynamics_application_variables.h"
#include "includes/cfd_variables.h"
#include "fluid_dynamics_application_variables.h"
#include "custom_utilities/potential_flow_utilities.h"

namespace Kratos
{
///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template <int Dim, int NumNodes>
Element::Pointer CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::Create(
    IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<CompressiblePerturbationPotentialFlowElement>(
        NewId, GetGeometry().Create(ThisNodes), pProperties);
    KRATOS_CATCH("");
}

template <int Dim, int NumNodes>
Element::Pointer CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::Create(
    IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<CompressiblePerturbationPotentialFlowElement>(NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

template <int Dim, int NumNodes>
Element::Pointer CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::Clone(
    IndexType NewId, NodesArrayType const& ThisNodes) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<CompressiblePerturbationPotentialFlowElement>(
        NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    KRATOS_CATCH("");
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    CalculateRightHandSide(rRightHandSideVector,rCurrentProcessInfo);
    CalculateLeftHandSide(rLeftHandSideMatrix,rCurrentProcessInfo);
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::CalculateRightHandSide(
    VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    const CompressiblePerturbationPotentialFlowElement& r_this = *this;
    const int wake = r_this.GetValue(WAKE);

    if (wake == 0) // Normal element (non-wake) - eventually an embedded
        CalculateRightHandSideNormalElement(rRightHandSideVector, rCurrentProcessInfo);
    else // Wake element
        CalculateRightHandSideWakeElement(rRightHandSideVector, rCurrentProcessInfo);
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    const CompressiblePerturbationPotentialFlowElement& r_this = *this;
    const int wake = r_this.GetValue(WAKE);

    if (wake == 0) // Normal element (non-wake) - eventually an embedded
        CalculateLeftHandSideNormalElement(rLeftHandSideMatrix, rCurrentProcessInfo);
    else // Wake element
        CalculateLeftHandSideWakeElement(rLeftHandSideMatrix, rCurrentProcessInfo);
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::EquationIdVector(
    EquationIdVectorType& rResult, const ProcessInfo& CurrentProcessInfo) const
{
    const CompressiblePerturbationPotentialFlowElement& r_this = *this;
    const int wake = r_this.GetValue(WAKE);

    if (wake == 0) // Normal element
    {
        if (rResult.size() != NumNodes)
            rResult.resize(NumNodes, false);

        const int kutta = r_this.GetValue(KUTTA);

        if (kutta == 0)
            GetEquationIdVectorNormalElement(rResult);
        else
            GetEquationIdVectorKuttaElement(rResult);
    }
    else // Wake element
    {
        if (rResult.size() != 2 * NumNodes)
            rResult.resize(2 * NumNodes, false);

        GetEquationIdVectorWakeElement(rResult);
    }
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::GetDofList(DofsVectorType& rElementalDofList,
                                                                 const ProcessInfo& CurrentProcessInfo) const
{
    const CompressiblePerturbationPotentialFlowElement& r_this = *this;
    const int wake = r_this.GetValue(WAKE);

    if (wake == 0) // Normal element
    {
        if (rElementalDofList.size() != NumNodes)
            rElementalDofList.resize(NumNodes);

        const int kutta = r_this.GetValue(KUTTA);

        if (kutta == 0)
            GetDofListNormalElement(rElementalDofList);
        else
            GetDofListKuttaElement(rElementalDofList);
    }
    else // wake element
    {
        if (rElementalDofList.size() != 2 * NumNodes)
            rElementalDofList.resize(2 * NumNodes);

        GetDofListWakeElement(rElementalDofList);
    }
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    bool active = true;
    if ((this)->IsDefined(ACTIVE))
        active = (this)->Is(ACTIVE);

    const CompressiblePerturbationPotentialFlowElement& r_this = *this;
    const int wake = r_this.GetValue(WAKE);

    if (wake != 0 && active == true)
    {
        ComputePotentialJump(rCurrentProcessInfo);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Inquiry

template <int Dim, int NumNodes>
int CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Generic geometry check
    int out = Element::Check(rCurrentProcessInfo);
    if (out != 0)
    {
        return out;
    }

    KRATOS_ERROR_IF(GetGeometry().Area() <= 0.0)
        << this->Id() << "Area cannot be less than or equal to 0" << std::endl;

    for (unsigned int i = 0; i < this->GetGeometry().size(); i++)
    {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY_POTENTIAL, this->GetGeometry()[i]);
    }

    return out;

    KRATOS_CATCH("");
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1)
        rValues.resize(1);

    if (rVariable == PRESSURE_COEFFICIENT)
    {
        rValues[0] = PotentialFlowUtilities::ComputePerturbationCompressiblePressureCoefficient<Dim, NumNodes>(*this, rCurrentProcessInfo);
    }
    else if (rVariable == DENSITY)
    {
        const array_1d<double, Dim> velocity = PotentialFlowUtilities::ComputePerturbedVelocity<Dim,NumNodes>(*this, rCurrentProcessInfo);
        const double local_mach_number_squared = PotentialFlowUtilities::ComputeLocalMachNumberSquared<Dim, NumNodes>(velocity, rCurrentProcessInfo);
        rValues[0] = PotentialFlowUtilities::ComputeDensity<Dim, NumNodes>(local_mach_number_squared, rCurrentProcessInfo);
    }
    else if (rVariable == MACH)
    {
        rValues[0] = PotentialFlowUtilities::ComputePerturbationLocalMachNumber<Dim, NumNodes>(*this, rCurrentProcessInfo);
    }
    else if (rVariable == SOUND_VELOCITY)
    {
        rValues[0] = PotentialFlowUtilities::ComputePerturbationLocalSpeedOfSound<Dim, NumNodes>(*this, rCurrentProcessInfo);
    }
    else if (rVariable == WAKE)
    {
        const CompressiblePerturbationPotentialFlowElement& r_this = *this;
        rValues[0] = r_this.GetValue(WAKE);
    }
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::CalculateOnIntegrationPoints(
    const Variable<int>& rVariable, std::vector<int>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1)
        rValues.resize(1);
    if (rVariable == TRAILING_EDGE)
        rValues[0] = this->GetValue(TRAILING_EDGE);
    else if (rVariable == KUTTA)
        rValues[0] = this->GetValue(KUTTA);
    else if (rVariable == WAKE)
        rValues[0] = this->GetValue(WAKE);
    else if (rVariable == ZERO_VELOCITY_CONDITION)
        rValues[0] = this->GetValue(ZERO_VELOCITY_CONDITION);
    else if (rVariable == TRAILING_EDGE_ELEMENT)
        rValues[0] = this->GetValue(TRAILING_EDGE_ELEMENT);
    else if (rVariable == DECOUPLED_TRAILING_EDGE_ELEMENT)
        rValues[0] = this->GetValue(DECOUPLED_TRAILING_EDGE_ELEMENT);
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable,
    std::vector<array_1d<double, 3>>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1)
        rValues.resize(1);
    if (rVariable == VELOCITY){
        const array_1d<double, 3> free_stream_velocity = rCurrentProcessInfo[FREE_STREAM_VELOCITY];
        array_1d<double, 3> v(3, 0.0);
        array_1d<double, Dim> vaux = PotentialFlowUtilities::ComputeVelocity<Dim, NumNodes>(*this);
        for (unsigned int k = 0; k < Dim; k++)
            v[k] = vaux[k] + free_stream_velocity[k];
        rValues[0] = v;
    }
    else if (rVariable == PERTURBATION_VELOCITY)
    {
        array_1d<double, 3> v(3, 0.0);
        array_1d<double, Dim> vaux = PotentialFlowUtilities::ComputeVelocity<Dim,NumNodes>(*this);
        for (unsigned int k = 0; k < Dim; k++)
            v[k] = vaux[k];
        rValues[0] = v;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Input and output

template <int Dim, int NumNodes>
std::string CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::Info() const
{
    std::stringstream buffer;
    buffer << "CompressiblePerturbationPotentialFlowElement #" << Id();
    return buffer.str();
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "CompressiblePerturbationPotentialFlowElement #" << Id();
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::PrintData(std::ostream& rOStream) const
{
    pGetGeometry()->PrintData(rOStream);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::GetWakeDistances(
    array_1d<double, NumNodes>& distances) const
{
    noalias(distances) = GetValue(WAKE_ELEMENTAL_DISTANCES);
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::GetEquationIdVectorNormalElement(
    EquationIdVectorType& rResult) const
{
    for (unsigned int i = 0; i < NumNodes; i++)
        rResult[i] = GetGeometry()[i].GetDof(VELOCITY_POTENTIAL).EquationId();
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::GetEquationIdVectorKuttaElement(
    EquationIdVectorType& rResult) const
{
    const auto& r_geometry = this->GetGeometry();
    // Kutta elements have only negative part
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        if (!r_geometry[i].GetValue(TRAILING_EDGE))
            rResult[i] = r_geometry[i].GetDof(VELOCITY_POTENTIAL).EquationId();
        else
            rResult[i] = r_geometry[i].GetDof(AUXILIARY_VELOCITY_POTENTIAL).EquationId();
    }
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::GetEquationIdVectorWakeElement(
    EquationIdVectorType& rResult) const
{
    array_1d<double, NumNodes> distances;
    GetWakeDistances(distances);

    // Positive part
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        if (distances[i] > 0.0)
            rResult[i] = GetGeometry()[i].GetDof(VELOCITY_POTENTIAL).EquationId();
        else
            rResult[i] =
                GetGeometry()[i].GetDof(AUXILIARY_VELOCITY_POTENTIAL, 0).EquationId();
    }

    // Negative part - sign is opposite to the previous case
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        if (distances[i] < 0.0)
            rResult[NumNodes + i] =
                GetGeometry()[i].GetDof(VELOCITY_POTENTIAL).EquationId();
        else
            rResult[NumNodes + i] =
                GetGeometry()[i].GetDof(AUXILIARY_VELOCITY_POTENTIAL).EquationId();
    }
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::GetDofListNormalElement(DofsVectorType& rElementalDofList) const
{
    for (unsigned int i = 0; i < NumNodes; i++)
        rElementalDofList[i] = GetGeometry()[i].pGetDof(VELOCITY_POTENTIAL);
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::GetDofListKuttaElement(DofsVectorType& rElementalDofList) const
{
    const auto& r_geometry = this->GetGeometry();
    // Kutta elements have only negative part
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        if (!r_geometry[i].GetValue(TRAILING_EDGE))
            rElementalDofList[i] = r_geometry[i].pGetDof(VELOCITY_POTENTIAL);
        else
            rElementalDofList[i] = r_geometry[i].pGetDof(AUXILIARY_VELOCITY_POTENTIAL);
    }
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::GetDofListWakeElement(DofsVectorType& rElementalDofList) const
{
    array_1d<double, NumNodes> distances;
    GetWakeDistances(distances);

    // Positive part
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        if (distances[i] > 0)
            rElementalDofList[i] = GetGeometry()[i].pGetDof(VELOCITY_POTENTIAL);
        else
            rElementalDofList[i] = GetGeometry()[i].pGetDof(AUXILIARY_VELOCITY_POTENTIAL);
    }

    // Negative part - sign is opposite to the previous case
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        if (distances[i] < 0)
            rElementalDofList[NumNodes + i] = GetGeometry()[i].pGetDof(VELOCITY_POTENTIAL);
        else
            rElementalDofList[NumNodes + i] =
                GetGeometry()[i].pGetDof(AUXILIARY_VELOCITY_POTENTIAL);
    }
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::CalculateLeftHandSideNormalElement(
    MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != NumNodes || rLeftHandSideMatrix.size2() != NumNodes)
        rLeftHandSideMatrix.resize(NumNodes, NumNodes, false);
    rLeftHandSideMatrix.clear();

    ElementalData<NumNodes, Dim> data;
    GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

    const array_1d<double, Dim> velocity = PotentialFlowUtilities::ComputePerturbedVelocity<Dim,NumNodes>(*this, rCurrentProcessInfo);

    const double local_mach_number_squared = PotentialFlowUtilities::ComputeLocalMachNumberSquared<Dim, NumNodes>(velocity, rCurrentProcessInfo);

    const double density = PotentialFlowUtilities::ComputeDensity<Dim, NumNodes>(local_mach_number_squared, rCurrentProcessInfo);

    const double DrhoDu2 = PotentialFlowUtilities::ComputeDensityDerivativeWRTVelocitySquared<Dim, NumNodes>(local_mach_number_squared, rCurrentProcessInfo);

    const BoundedVector<double, NumNodes> DNV = prod(data.DN_DX, velocity);

    noalias(rLeftHandSideMatrix) += data.vol * density * prod(data.DN_DX, trans(data.DN_DX));

    const double max_velocity_squared = PotentialFlowUtilities::ComputeMaximumVelocitySquared<Dim, NumNodes>(rCurrentProcessInfo);

    const double local_velocity_squared = inner_prod(velocity, velocity);

    if (local_velocity_squared < max_velocity_squared){
        noalias(rLeftHandSideMatrix) += data.vol * 2 * DrhoDu2 * outer_prod(DNV, trans(DNV));
    }
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::CalculateRightHandSideNormalElement(
    VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    if (rRightHandSideVector.size() != NumNodes)
        rRightHandSideVector.resize(NumNodes, false);
    rRightHandSideVector.clear();

    ElementalData<NumNodes, Dim> data;

    // Calculate shape functions
    GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

    const CompressiblePerturbationPotentialFlowElement& r_this = *this;
    array_1d<double, Dim> velocity = PotentialFlowUtilities::ComputePerturbedVelocity<Dim,NumNodes>(r_this, rCurrentProcessInfo);

    const double max_velocity_squared = PotentialFlowUtilities::ComputeMaximumVelocitySquared<Dim, NumNodes>(rCurrentProcessInfo);

    if (this->Id()==4){
        KRATOS_WATCH(std::sqrt(max_velocity_squared))
    }

    // const double local_velocity_squared = inner_prod(velocity, velocity);

    // if (local_velocity_squared > max_velocity_squared){
    //     velocity *= std::sqrt(max_velocity_squared / local_velocity_squared);
    // }

    const double local_mach_number_squared = PotentialFlowUtilities::ComputeLocalMachNumberSquared<Dim, NumNodes>(velocity, rCurrentProcessInfo);
    const double density = PotentialFlowUtilities::ComputeDensity<Dim, NumNodes>(local_mach_number_squared, rCurrentProcessInfo);

    noalias(rRightHandSideVector) = - data.vol * density * prod(data.DN_DX, velocity);
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::CalculateLeftHandSideWakeElement(
    MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    // Note that the lhs and rhs have double the size
    if (rLeftHandSideMatrix.size1() != 2 * NumNodes ||
        rLeftHandSideMatrix.size2() != 2 * NumNodes)
        rLeftHandSideMatrix.resize(2 * NumNodes, 2 * NumNodes, false);
    rLeftHandSideMatrix.clear();

    ElementalData<NumNodes, Dim> data;

    // Calculate shape functions
    GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);
    GetWakeDistances(data.distances);

    const array_1d<double, TDim> upper_velocity = PotentialFlowUtilities::ComputePerturbedVelocity<Dim,NumNodes>(*this, rCurrentProcessInfo);

    const double local_mach_number_squared = PotentialFlowUtilities::ComputeLocalMachNumberSquared<Dim, NumNodes>(upper_velocity, rCurrentProcessInfo);

    const double density = PotentialFlowUtilities::ComputeDensity<Dim, NumNodes>(local_mach_number_squared, rCurrentProcessInfo);

    const double DrhoDu2 = PotentialFlowUtilities::ComputeDensityDerivativeWRTVelocitySquared<Dim, NumNodes>(local_mach_number_squared, rCurrentProcessInfo);

    const BoundedVector<double, NumNodes> DNV = prod(data.DN_DX, upper_velocity);

    const double max_velocity_squared = PotentialFlowUtilities::ComputeMaximumVelocitySquared<Dim, NumNodes>(rCurrentProcessInfo);

    const double local_velocity_squared = inner_prod(upper_velocity, upper_velocity);

    BoundedMatrix<double, NumNodes, NumNodes> lhs_total = data.vol * density * prod(data.DN_DX, trans(data.DN_DX));

    if (local_velocity_squared < max_velocity_squared){
        lhs_total += data.vol * 2 * DrhoDu2 * outer_prod(DNV, trans(DNV));
    }

    const array_1d<double, 3> free_stream_velocity = rCurrentProcessInfo[FREE_STREAM_VELOCITY];
    array_1d<double, Dim> lower_velocity = PotentialFlowUtilities::ComputeVelocityLowerWakeElement<Dim,NumNodes>(*this);

    for (unsigned int i = 0; i < Dim; i++){
        lower_velocity[i] += free_stream_velocity[i];
    }

    const double local_mach_number_squared_lower = PotentialFlowUtilities::ComputeLocalMachNumberSquared<Dim, NumNodes>(lower_velocity, rCurrentProcessInfo);
    const double density_lower = PotentialFlowUtilities::ComputeDensity<Dim, NumNodes>(local_mach_number_squared_lower, rCurrentProcessInfo);
    const double DrhoDu2_lower = PotentialFlowUtilities::ComputeDensityDerivativeWRTVelocitySquared<Dim, NumNodes>(local_mach_number_squared_lower, rCurrentProcessInfo);

    const BoundedVector<double, NumNodes> DNV_lower = prod(data.DN_DX, lower_velocity);

    BoundedMatrix<double, NumNodes, NumNodes> lhs_total_lower = data.vol * density_lower * prod(data.DN_DX, trans(data.DN_DX));

    const double local_velocity_squared_lower = inner_prod(lower_velocity, lower_velocity);

    if (local_velocity_squared_lower < max_velocity_squared){
        lhs_total_lower += data.vol * 2 * DrhoDu2_lower * outer_prod(DNV_lower, trans(DNV_lower));
    }

    // Wake condition matrix
    BoundedMatrix<double, Dim, Dim> condition_matrix = IdentityMatrix(Dim,Dim);
    condition_matrix(0,0) = 1.0;
    condition_matrix(1,1) = 0.0;
    condition_matrix(2,2) = 0.0;

    auto xzfilter = prod(condition_matrix, trans(data.DN_DX));

    BoundedMatrix<double, NumNodes, NumNodes> lhs_wake_condition = ZeroMatrix(NumNodes, NumNodes);

    const double free_stream_density = rCurrentProcessInfo[FREE_STREAM_DENSITY];
    for (unsigned int i = 0; i < NumNodes; i++){
        for(unsigned int j = 0; j < NumNodes; j++){
            for(unsigned int k = 0; k < Dim; k++){
                lhs_wake_condition(i,j) += data.vol*free_stream_density*data.DN_DX(i,k)*xzfilter(k,j);
            }
        }
    }

    const double velocity_upper_2 = inner_prod(upper_velocity, upper_velocity);
    const double velocity_lower_2 = inner_prod(lower_velocity, lower_velocity);
    const double diff_upper = velocity_upper_2 - velocity_lower_2;
    const double diff_lower = velocity_lower_2 - velocity_upper_2;

    const double free_stream_mach = rCurrentProcessInfo[FREE_STREAM_MACH];
    const double free_stream_speed_of_sound = rCurrentProcessInfo[SOUND_VELOCITY];
    const double u_inf = free_stream_mach * free_stream_speed_of_sound;
    const double stabilization_upper = 0.5 * std::pow(data.vol, 1.0/3.0) * inner_prod(free_stream_velocity, upper_velocity) / u_inf;
    const double stabilization_lower = 0.5 * std::pow(data.vol, 1.0/3.0) * inner_prod(free_stream_velocity, lower_velocity) / u_inf;

    BoundedVector<double, NumNodes> dU2dPhi_upper = 2*data.vol*prod(data.DN_DX, upper_velocity);
    BoundedVector<double, NumNodes> dU2dPhi_lower = 2*data.vol*prod(data.DN_DX, lower_velocity);

    BoundedMatrix<double, NumNodes, NumNodes> lhs_pressure_upper_up = ZeroMatrix(NumNodes, NumNodes);
    BoundedMatrix<double, NumNodes, NumNodes> lhs_pressure_upper_low = ZeroMatrix(NumNodes, NumNodes);
    BoundedMatrix<double, NumNodes, NumNodes> lhs_pressure_lower_low = ZeroMatrix(NumNodes, NumNodes);
    BoundedMatrix<double, NumNodes, NumNodes> lhs_pressure_lower_up = ZeroMatrix(NumNodes, NumNodes);

    const BoundedVector<double, NumNodes> lhs_stabilization_upper =  0.5 * std::pow(data.vol, 1.0/3.0) * diff_upper * prod(data.DN_DX, free_stream_velocity) / u_inf * data.vol;
    const BoundedVector<double, NumNodes> lhs_stabilization_lower =  0.5 * std::pow(data.vol, 1.0/3.0) * diff_lower * prod(data.DN_DX, free_stream_velocity) / u_inf * data.vol;

    for(unsigned int row = 0; row < NumNodes; ++row){
        for(unsigned int column = 0; column < NumNodes; ++column){
            lhs_pressure_upper_up(row, column) += lhs_stabilization_upper[column];
            lhs_pressure_upper_up(row, column) += (data.N[row] + stabilization_upper) * dU2dPhi_upper[column];

            // lhs_pressure_upper_low(row, column) += lhs_stabilization_upper[column];
            lhs_pressure_upper_low(row, column) += (data.N[row] + stabilization_upper) * dU2dPhi_lower[column];

            lhs_pressure_lower_low(row, column) += lhs_stabilization_lower[column];
            lhs_pressure_lower_low(row, column) += (data.N[row] + stabilization_lower) * dU2dPhi_lower[column];

            // lhs_pressure_lower_up(row, column) += lhs_stabilization_lower[column];
            lhs_pressure_lower_up(row, column) += (data.N[row] + stabilization_lower) * dU2dPhi_upper[column];
        }
    }

    BoundedMatrix<double, NumNodes, NumNodes> lhs_positive = ZeroMatrix(NumNodes, NumNodes);
    BoundedMatrix<double, NumNodes, NumNodes> lhs_negative = ZeroMatrix(NumNodes, NumNodes);

    CalculateLeftHandSideSubdividedElement(lhs_positive, lhs_negative, rCurrentProcessInfo);

    if (this->Is(STRUCTURE)){

        const auto& r_geometry = this->GetGeometry();
        for (unsigned int row = 0; row < NumNodes; ++row){
            // The TE node takes the contribution of the subdivided element and
            // we do not apply the wake condition on the TE node
            if (r_geometry[row].GetValue(TRAILING_EDGE)){
                for (unsigned int column = 0; column < NumNodes; ++column){
                    // Conservation of mass
                    rLeftHandSideMatrix(row, column) = lhs_positive(row, column);
                    rLeftHandSideMatrix(row + NumNodes, column + NumNodes) = lhs_negative(row, column);

                    // // Wake condition below
                    // rLeftHandSideMatrix(row + NumNodes, column + NumNodes) = lhs_wake_condition(row, column); // Diagonal
                    // rLeftHandSideMatrix(row + NumNodes, column) = -lhs_wake_condition(row, column); // Off diagonal
                }
            }
            else{
                // Applying wake condition on the AUXILIARY_VELOCITY_POTENTIAL dofs
                if (data.distances[row] < 0.0){
                    for (unsigned int column = 0; column < NumNodes; ++column){
                        // Conservation of mass
                        rLeftHandSideMatrix(row + NumNodes, column + NumNodes) = lhs_negative(row, column);
                        // Wake condition
                        rLeftHandSideMatrix(row, column) = lhs_wake_condition(row, column); // Diagonal
                        rLeftHandSideMatrix(row, column + NumNodes) = -lhs_wake_condition(row, column); // Off diagonal
                        // rLeftHandSideMatrix(row, column) = lhs_pressure_upper_up(row, column); // Diagonal
                        // rLeftHandSideMatrix(row, column + NumNodes) = -lhs_pressure_upper_low(row, column); // Off diagonal
                    }
                }
                else{ // else if (data.distances[row] > 0.0)
                    for (unsigned int column = 0; column < NumNodes; ++column){
                        // Conservation of mass
                        rLeftHandSideMatrix(row, column) = lhs_positive(row, column);
                        // Wake condition
                        rLeftHandSideMatrix(row + NumNodes, column + NumNodes) = lhs_wake_condition(row, column); // Diagonal
                        rLeftHandSideMatrix(row + NumNodes, column) = -lhs_wake_condition(row, column); // Off diagonal
                        // rLeftHandSideMatrix(row + NumNodes, column + NumNodes) = lhs_pressure_lower_low(row, column); // Diagonal
                        // rLeftHandSideMatrix(row + NumNodes, column) = -lhs_pressure_lower_up(row, column); // Off diagonal
                    }
                }
            }
        }
    }
    else{
        for (unsigned int row = 0; row < NumNodes; ++row){
            // Applying wake condition on the AUXILIARY_VELOCITY_POTENTIAL dofs
            if (data.distances[row] < 0.0){
                for (unsigned int column = 0; column < NumNodes; ++column){
                    // Conservation of mass
                    rLeftHandSideMatrix(row + NumNodes, column + NumNodes) = lhs_negative(row, column);
                    // Wake condition
                    rLeftHandSideMatrix(row, column) = lhs_wake_condition(row, column); // Diagonal
                    rLeftHandSideMatrix(row, column + NumNodes) = -lhs_wake_condition(row, column); // Off diagonal
                    // rLeftHandSideMatrix(row, column) = lhs_pressure_upper_up(row, column); // Diagonal
                    // rLeftHandSideMatrix(row, column + NumNodes) = -lhs_pressure_upper_low(row, column); // Off diagonal
                }
            }
            else{ // else if (data.distances[row] > 0.0)
                for (unsigned int column = 0; column < NumNodes; ++column){
                    // Conservation of mass
                    rLeftHandSideMatrix(row, column) = lhs_positive(row, column);
                    // Wake condition
                    rLeftHandSideMatrix(row + NumNodes, column + NumNodes) = lhs_wake_condition(row, column); // Diagonal
                    rLeftHandSideMatrix(row + NumNodes, column) = -lhs_wake_condition(row, column); // Off diagonal
                    //  rLeftHandSideMatrix(row + NumNodes, column + NumNodes) = lhs_pressure_lower_low(row, column); // Diagonal
                    // rLeftHandSideMatrix(row + NumNodes, column) = -lhs_pressure_lower_up(row, column); // Off diagonal
                }
            }
        }
    }
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::CalculateRightHandSideWakeElement(
    VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    // Note that the rhs has double the size
    if (rRightHandSideVector.size() != 2 * NumNodes)
        rRightHandSideVector.resize(2 * NumNodes, false);
    rRightHandSideVector.clear();

    ElementalData<NumNodes, Dim> data;

    // Calculate shape functions
    const auto& r_geometry = this->GetGeometry();
    GeometryUtils::CalculateGeometryData(r_geometry, data.DN_DX, data.N, data.vol);
    GetWakeDistances(data.distances);

    const array_1d<double, 3> free_stream_velocity = rCurrentProcessInfo[FREE_STREAM_VELOCITY];
    array_1d<double, Dim> upper_velocity = PotentialFlowUtilities::ComputeVelocityUpperWakeElement<Dim,NumNodes>(*this);
    array_1d<double, Dim> lower_velocity = PotentialFlowUtilities::ComputeVelocityLowerWakeElement<Dim,NumNodes>(*this);

    for (unsigned int i = 0; i < Dim; i++){
        upper_velocity[i] += free_stream_velocity[i];
        lower_velocity[i] += free_stream_velocity[i];
    }

    // const double max_velocity_squared = PotentialFlowUtilities::ComputeMaximumVelocitySquared<Dim, NumNodes>(rCurrentProcessInfo);

    // const double upper_local_velocity_squared = inner_prod(upper_velocity, upper_velocity);
    // const double lower_local_velocity_squared = inner_prod(lower_velocity, lower_velocity);

    // if (upper_local_velocity_squared > max_velocity_squared || lower_local_velocity_squared > max_velocity_squared){
    //     upper_velocity *= std::sqrt(max_velocity_squared / upper_local_velocity_squared);
    //     lower_velocity *= std::sqrt(max_velocity_squared / lower_local_velocity_squared);
    // }

    const double local_mach_number_squared_upper = PotentialFlowUtilities::ComputeLocalMachNumberSquared<Dim, NumNodes>(upper_velocity, rCurrentProcessInfo);
    const double density_upper = PotentialFlowUtilities::ComputeDensity<Dim, NumNodes>(local_mach_number_squared_upper, rCurrentProcessInfo);

    const double local_mach_number_squared_lower = PotentialFlowUtilities::ComputeLocalMachNumberSquared<Dim, NumNodes>(lower_velocity, rCurrentProcessInfo);
    const double density_lower = PotentialFlowUtilities::ComputeDensity<Dim, NumNodes>(local_mach_number_squared_lower, rCurrentProcessInfo);

    const array_1d<double, Dim> diff_velocity = upper_velocity - lower_velocity;

    const BoundedVector<double, NumNodes> upper_rhs = - data.vol * density_upper * prod(data.DN_DX, upper_velocity);
    const BoundedVector<double, NumNodes> lower_rhs = - data.vol * density_lower * prod(data.DN_DX, lower_velocity);

    // Wake condition matrix
    BoundedMatrix<double, Dim, Dim> condition_matrix = IdentityMatrix(Dim,Dim);
    condition_matrix(0,0) = 1.0;
    condition_matrix(1,1) = 0.0;
    condition_matrix(2,2) = 0.0;

    const auto xzfilter = prod(data.DN_DX, condition_matrix);

    const double free_stream_density = rCurrentProcessInfo[FREE_STREAM_DENSITY];
    const BoundedVector<double, NumNodes> wake_rhs_upper = - data.vol * free_stream_density * prod(xzfilter, diff_velocity);
    const BoundedVector<double, NumNodes> wake_rhs_lower = data.vol * free_stream_density * prod(xzfilter, diff_velocity);

    // const double velocity_upper_2 = inner_prod(upper_velocity, upper_velocity);
    // const double velocity_lower_2 = inner_prod(lower_velocity, lower_velocity);


    // const double free_stream_mach = rCurrentProcessInfo[FREE_STREAM_MACH];
    // const double free_stream_speed_of_sound = rCurrentProcessInfo[SOUND_VELOCITY];
    // const double u_inf = free_stream_mach * free_stream_speed_of_sound;
    // const double stabilization_upper = 0.5 * std::pow(data.vol, 1.0/3.0) * inner_prod(free_stream_velocity, upper_velocity) / u_inf;
    // const double stabilization_lower = 0.5 * std::pow(data.vol, 1.0/3.0) * inner_prod(free_stream_velocity, lower_velocity) / u_inf;

    // BoundedVector<double, NumNodes> wake_rhs_upper = ZeroVector(NumNodes);
    // BoundedVector<double, NumNodes> wake_rhs_lower = ZeroVector(NumNodes);
    // for(unsigned int row = 0; row < NumNodes; ++row){
    //     wake_rhs_upper(row) = - (data.N[row] + stabilization_upper) * (velocity_upper_2 - velocity_lower_2) * data.vol;
    //     wake_rhs_lower(row) = - (data.N[row] + stabilization_lower) * (velocity_lower_2 - velocity_upper_2) * data.vol;
    // }

    double upper_vol = 0.0;
    double lower_vol = 0.0;

    CalculateVolumesSubdividedElement(upper_vol, lower_vol, rCurrentProcessInfo);

    if (this->Is(STRUCTURE)){

        for (unsigned int rRow = 0; rRow < NumNodes; ++rRow){
            if (r_geometry[rRow].GetValue(TRAILING_EDGE)){
                rRightHandSideVector[rRow] = upper_rhs(rRow)*upper_vol/data.vol;
                rRightHandSideVector[rRow + NumNodes] = lower_rhs(rRow)*lower_vol/data.vol;
            }
            else{
                if (data.distances[rRow] > 0.0){
                    rRightHandSideVector[rRow] = upper_rhs(rRow)*upper_vol/data.vol;
                    rRightHandSideVector[rRow + TNumNodes] = wake_rhs_lower(rRow);
                }
                else{
                    rRightHandSideVector[rRow] = wake_rhs_upper(rRow);
                    rRightHandSideVector[rRow + NumNodes] = lower_rhs(rRow)*lower_vol/data.vol;
                }
            }
        }
    }
    else{
        for (unsigned int rRow = 0; rRow < NumNodes; ++rRow){
            if (data.distances[rRow] > 0.0){
                rRightHandSideVector[rRow] = upper_rhs(rRow)*upper_vol/data.vol;
                rRightHandSideVector[rRow + TNumNodes] = wake_rhs_lower(rRow);
            }
            else{
                rRightHandSideVector[rRow] = wake_rhs_upper(rRow);
                rRightHandSideVector[rRow + NumNodes] = lower_rhs(rRow)*lower_vol/data.vol;
            }
        }
    }
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::CalculateLeftHandSideSubdividedElement(
    BoundedMatrix<double, NumNodes, NumNodes>& lhs_positive,
    BoundedMatrix<double, NumNodes, NumNodes>& lhs_negative,
    const ProcessInfo& rCurrentProcessInfo)
{
    ElementalData<NumNodes, Dim> data;

    // Calculate shape functions
    GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

    GetWakeDistances(data.distances);

    // // Subdivide the element
    // constexpr unsigned int nvolumes = 3 * (Dim - 1);
    // BoundedMatrix<double, NumNodes, Dim> Points;
    // array_1d<double, nvolumes> PartitionsSign;
    // BoundedMatrix<double, nvolumes, NumNodes> GPShapeFunctionValues;
    // array_1d<double, nvolumes> Volumes;
    // std::vector<Matrix> GradientsValue(nvolumes);
    // BoundedMatrix<double, nvolumes, 2> NEnriched;
    // for (unsigned int i = 0; i < GradientsValue.size(); ++i)
    //     GradientsValue[i].resize(2, Dim, false);
    // for (unsigned int i = 0; i < NumNodes; ++i)
    // {
    //     const array_1d<double, 3>& coords = GetGeometry()[i].Coordinates();
    //     for (unsigned int k = 0; k < Dim; ++k)
    //     {
    //         Points(i, k) = coords[k];
    //     }
    // }

    // const unsigned int nsubdivisions = EnrichmentUtilities::CalculateEnrichedShapeFuncions(
    //     Points, data.DN_DX, data.distances, Volumes, GPShapeFunctionValues,
    //     PartitionsSign, GradientsValue, NEnriched);

    ModifiedShapeFunctions::Pointer p_calculator =
        EmbeddedDiscontinuousInternals::GetShapeFunctionCalculator<Dim, NumNodes>(
            *this,
            data.distances);

    Matrix positive_side_sh_func;
    ModifiedShapeFunctions::ShapeFunctionsGradientsType positive_side_sh_func_gradients;
    Vector positive_side_weights;

    p_calculator->ComputePositiveSideShapeFunctionsAndGradientsValues(
        positive_side_sh_func, positive_side_sh_func_gradients,
        positive_side_weights, GeometryData::GI_GAUSS_1);

    Matrix negative_side_sh_func;
    ModifiedShapeFunctions::ShapeFunctionsGradientsType negative_side_sh_func_gradients;
    Vector negative_side_weights;

    p_calculator->ComputeNegativeSideShapeFunctionsAndGradientsValues(
        negative_side_sh_func, negative_side_sh_func_gradients,
        negative_side_weights, GeometryData::GI_GAUSS_1);

    // Upper part
    const array_1d<double, TDim> velocity = PotentialFlowUtilities::ComputePerturbedVelocity<Dim,NumNodes>(*this, rCurrentProcessInfo);

    const double max_velocity_squared = PotentialFlowUtilities::ComputeMaximumVelocitySquared<Dim, NumNodes>(rCurrentProcessInfo);

    const double local_velocity_squared = inner_prod(velocity, velocity);

    const double local_mach_number_squared = PotentialFlowUtilities::ComputeLocalMachNumberSquared<Dim, NumNodes>(velocity, rCurrentProcessInfo);

    const double density = PotentialFlowUtilities::ComputeDensity<Dim, NumNodes>(local_mach_number_squared, rCurrentProcessInfo);

    const double DrhoDu2 = PotentialFlowUtilities::ComputeDensityDerivativeWRTVelocitySquared<Dim, NumNodes>(local_mach_number_squared, rCurrentProcessInfo);

    // const BoundedVector<double, NumNodes> DNV = prod(data.DN_DX, velocity);

    // Lower part
    const array_1d<double, 3> free_stream_velocity = rCurrentProcessInfo[FREE_STREAM_VELOCITY];
    array_1d<double, Dim> lower_velocity = PotentialFlowUtilities::ComputeVelocityLowerWakeElement<Dim,NumNodes>(*this);

    for (unsigned int i = 0; i < Dim; i++){
        lower_velocity[i] += free_stream_velocity[i];
    }

    const double local_mach_number_squared_lower = PotentialFlowUtilities::ComputeLocalMachNumberSquared<Dim, NumNodes>(lower_velocity, rCurrentProcessInfo);
    const double density_lower = PotentialFlowUtilities::ComputeDensity<Dim, NumNodes>(local_mach_number_squared_lower, rCurrentProcessInfo);
    const double DrhoDu2_lower = PotentialFlowUtilities::ComputeDensityDerivativeWRTVelocitySquared<Dim, NumNodes>(local_mach_number_squared_lower, rCurrentProcessInfo);

    // const BoundedVector<double, NumNodes> DNV_lower = prod(data.DN_DX, lower_velocity);

    const double local_velocity_squared_lower = inner_prod(lower_velocity, lower_velocity);

    for (unsigned int i_gauss = 0;
         i_gauss < positive_side_sh_func_gradients.size(); i_gauss++) {
        const BoundedMatrix<double, NumNodes, Dim> DN_DX = positive_side_sh_func_gradients(i_gauss);
        const BoundedVector<double, NumNodes> DNV = prod(DN_DX, velocity);
        noalias(lhs_positive) += positive_side_weights(i_gauss) * density * prod(DN_DX, trans(DN_DX));
        if(local_velocity_squared < max_velocity_squared){
            noalias(lhs_positive) += positive_side_weights(i_gauss) * 2 * DrhoDu2 * outer_prod(DNV, trans(DNV));
        }
    }

    for (unsigned int i_gauss = 0;
         i_gauss < negative_side_sh_func_gradients.size(); i_gauss++) {
        const BoundedMatrix<double, NumNodes, Dim> DN_DX = negative_side_sh_func_gradients(i_gauss);
        const BoundedVector<double, NumNodes> DNV_lower = prod(DN_DX, lower_velocity);
        noalias(lhs_negative) += negative_side_weights(i_gauss) * density_lower * prod(DN_DX, trans(DN_DX));
        if(local_velocity_squared_lower < max_velocity_squared){
            noalias(lhs_negative) += negative_side_weights(i_gauss) * 2 * DrhoDu2_lower * outer_prod(DNV_lower, trans(DNV_lower));
        }
    }

    // // Compute the lhs and rhs that would correspond to it being divided
    // for (unsigned int i = 0; i < nsubdivisions; ++i)
    // {
    //     if (PartitionsSign[i] > 0)
    //     {
    //         noalias(lhs_positive) += Volumes[i] * density * prod(data.DN_DX, trans(data.DN_DX));
    //         if(local_velocity_squared < max_velocity_squared){
    //             noalias(lhs_positive) += Volumes[i] * 2 * DrhoDu2 * outer_prod(DNV, trans(DNV));
    //         }
    //     }
    //     else
    //     {
    //         noalias(lhs_negative) += Volumes[i] * density_lower * prod(data.DN_DX, trans(data.DN_DX));
    //         if(local_velocity_squared_lower < max_velocity_squared){
    //             noalias(lhs_negative) += Volumes[i] * 2 * DrhoDu2_lower * outer_prod(DNV_lower, trans(DNV_lower));
    //         }
    //     }
    // }
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::CalculateVolumesSubdividedElement(
    double& rUpper_vol,
    double& rLower_vol,
    const ProcessInfo& rCurrentProcessInfo)
{
    ElementalData<NumNodes, Dim> data;

    // Calculate shape functions
    GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

    GetWakeDistances(data.distances);

    // // Subdivide the element
    // constexpr unsigned int nvolumes = 3 * (Dim - 1);
    // BoundedMatrix<double, NumNodes, Dim> Points;
    // array_1d<double, nvolumes> PartitionsSign;
    // BoundedMatrix<double, nvolumes, NumNodes> GPShapeFunctionValues;
    // array_1d<double, nvolumes> Volumes;
    // std::vector<Matrix> GradientsValue(nvolumes);
    // BoundedMatrix<double, nvolumes, 2> NEnriched;
    // for (unsigned int i = 0; i < GradientsValue.size(); ++i)
    //     GradientsValue[i].resize(2, Dim, false);
    // for (unsigned int i = 0; i < NumNodes; ++i)
    // {
    //     const array_1d<double, 3>& coords = GetGeometry()[i].Coordinates();
    //     for (unsigned int k = 0; k < Dim; ++k)
    //     {
    //         Points(i, k) = coords[k];
    //     }
    // }

    // const unsigned int nsubdivisions = EnrichmentUtilities::CalculateEnrichedShapeFuncions(
    //     Points, data.DN_DX, data.distances, Volumes, GPShapeFunctionValues,
    //     PartitionsSign, GradientsValue, NEnriched);

    ModifiedShapeFunctions::Pointer p_calculator =
        EmbeddedDiscontinuousInternals::GetShapeFunctionCalculator<Dim, NumNodes>(
            *this,
            data.distances);

    Matrix positive_side_sh_func;
    ModifiedShapeFunctions::ShapeFunctionsGradientsType positive_side_sh_func_gradients;
    Vector positive_side_weights;

    p_calculator->ComputePositiveSideShapeFunctionsAndGradientsValues(
        positive_side_sh_func, positive_side_sh_func_gradients,
        positive_side_weights, GeometryData::GI_GAUSS_1);

    Matrix negative_side_sh_func;
    ModifiedShapeFunctions::ShapeFunctionsGradientsType negative_side_sh_func_gradients;
    Vector negative_side_weights;

    p_calculator->ComputeNegativeSideShapeFunctionsAndGradientsValues(
        negative_side_sh_func, negative_side_sh_func_gradients,
        negative_side_weights, GeometryData::GI_GAUSS_1);

    for (unsigned int i = 0; i < positive_side_sh_func_gradients.size(); ++i){
        rUpper_vol += positive_side_weights(i);
    }

    for (unsigned int i = 0; i < negative_side_sh_func_gradients.size(); ++i){
        rLower_vol += negative_side_weights(i);
    }

    // // Compute the volumes that would correspond to it being divided
    // for (unsigned int i = 0; i < nsubdivisions; ++i)
    // {
    //     if (PartitionsSign[i] > 0){
    //         rUpper_vol += Volumes[i];
    //     }
    //     else{
    //         rLower_vol += Volumes[i];
    //     }
    // }
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::ComputeLHSGaussPointContribution(
    const double weight, BoundedMatrix<double, NumNodes, NumNodes>& lhs, const ElementalData<NumNodes, Dim>& data) const
{
    noalias(lhs) += weight * prod(data.DN_DX, trans(data.DN_DX));
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::AssignLeftHandSideSubdividedElement(
    Matrix& rLeftHandSideMatrix,
    BoundedMatrix<double, NumNodes, NumNodes>& lhs_positive,
    BoundedMatrix<double, NumNodes, NumNodes>& lhs_negative,
    const BoundedMatrix<double, NumNodes, NumNodes>& lhs_total,
    const ElementalData<NumNodes, Dim>& data) const
{
    const auto& r_geometry = this->GetGeometry();
    for (unsigned int i = 0; i < NumNodes; ++i)
    {
        // The TE node takes the contribution of the subdivided element and
        // we do not apply the wake condition on the TE node
        if (r_geometry[i].GetValue(TRAILING_EDGE)){
            for (unsigned int j = 0; j < NumNodes; ++j){
                rLeftHandSideMatrix(i, j) = lhs_positive(i, j);
                rLeftHandSideMatrix(i + NumNodes, j + NumNodes) = lhs_negative(i, j);
            }
        }
        else{
            AssignLeftHandSideWakeNode(rLeftHandSideMatrix, lhs_total, data, i);
        }
    }
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::AssignLeftHandSideWakeElement(
    MatrixType& rLeftHandSideMatrix,
    const BoundedMatrix<double, NumNodes, NumNodes>& lhs_total,
    const ElementalData<NumNodes, Dim>& data) const
{
    for (unsigned int row = 0; row < NumNodes; ++row)
        AssignLeftHandSideWakeNode(rLeftHandSideMatrix, lhs_total, data, row);
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::AssignLeftHandSideWakeNode(
    MatrixType& rLeftHandSideMatrix,
    const BoundedMatrix<double, NumNodes, NumNodes>& lhs_total,
    const ElementalData<NumNodes, Dim>& data,
    unsigned int& row) const
{
    // Filling the diagonal blocks (i.e. decoupling upper and lower dofs)
    for (unsigned int column = 0; column < NumNodes; ++column)
    {
        rLeftHandSideMatrix(row, column) = lhs_total(row, column);
        rLeftHandSideMatrix(row + NumNodes, column + NumNodes) = lhs_total(row, column);
    }

    // Applying wake condition on the AUXILIARY_VELOCITY_POTENTIAL dofs
    if (data.distances[row] < 0.0)
        for (unsigned int column = 0; column < NumNodes; ++column)
            rLeftHandSideMatrix(row, column + NumNodes) = -lhs_total(row, column); // Side 1
    else if (data.distances[row] > 0.0)
        for (unsigned int column = 0; column < NumNodes; ++column)
            rLeftHandSideMatrix(row + NumNodes, column) = -lhs_total(row, column); // Side 2
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::AssignRightHandSideWakeNode(
    VectorType& rRightHandSideVector,
    const BoundedVector<double, NumNodes>& rUpper_rhs,
    const BoundedVector<double, NumNodes>& rLower_rhs,
    const BoundedVector<double, NumNodes>& rWake_rhs,
    const ElementalData<NumNodes, Dim>& rData,
    unsigned int& rRow) const
{
    if (rData.distances[rRow] > 0.0){
        rRightHandSideVector[rRow] = rUpper_rhs(rRow);
        rRightHandSideVector[rRow + TNumNodes] = -rWake_rhs(rRow);
    }
    else{
        rRightHandSideVector[rRow] = rWake_rhs(rRow);
        rRightHandSideVector[rRow + NumNodes] = rLower_rhs(rRow);
    }
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::ComputePotentialJump(const ProcessInfo& rCurrentProcessInfo)
{
    const array_1d<double, 3>& vinfinity = rCurrentProcessInfo[FREE_STREAM_VELOCITY];
    const double vinfinity_norm = sqrt(inner_prod(vinfinity, vinfinity));

    array_1d<double, NumNodes> distances;
    GetWakeDistances(distances);

    auto r_geometry = GetGeometry();
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        double aux_potential = r_geometry[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL);
        double potential = r_geometry[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL);
        double potential_jump = aux_potential - potential;

        if (distances[i] > 0)
        {
            r_geometry[i].SetLock();
            r_geometry[i].SetValue(POTENTIAL_JUMP, -2.0 / vinfinity_norm * (potential_jump));
            r_geometry[i].UnSetLock();
        }
        else
        {
            r_geometry[i].SetLock();
            r_geometry[i].SetValue(POTENTIAL_JUMP, 2.0 / vinfinity_norm * (potential_jump));
            r_geometry[i].UnSetLock();
        }
    }
}

// serializer

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Helper functions for template specialization
///////////////////////////////////////////////////////////////////////////////////////////////////

namespace EmbeddedDiscontinuousInternals {

template <>
ModifiedShapeFunctions::Pointer GetShapeFunctionCalculator<3, 4>(const Element& rElement, const Vector& rElementalDistances)
{
    return Kratos::make_shared<Tetrahedra3D4ModifiedShapeFunctions>(rElement.pGetGeometry(), rElementalDistances);
}

}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Template class instantiation

template class CompressiblePerturbationPotentialFlowElement<2, 3>;
template class CompressiblePerturbationPotentialFlowElement<3, 4>;

} // namespace Kratos
