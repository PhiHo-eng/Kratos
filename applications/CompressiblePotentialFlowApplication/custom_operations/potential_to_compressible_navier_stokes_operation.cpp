//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Marco Antonio Zuñiga Perez
//
//

// System includes

// External includes

// Project includes
#include "utilities/parallel_utilities.h"
#include "potential_to_compressible_navier_stokes_operation.h"
#include "fluid_dynamics_application_variables.h"
#include "compressible_potential_flow_application_variables.h"
#include "custom_processes/compute_nodal_value_process.h"

namespace Kratos
{

PotentialToCompressibleNavierStokesOperation::PotentialToCompressibleNavierStokesOperation(
    Model& rModel,
    Parameters ModelParameters)
        : mpModel(&rModel)
        , mParameters(ModelParameters)
{
    mParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());
}

Operation::Pointer PotentialToCompressibleNavierStokesOperation::Create(Model &rModel,
    Parameters Parameters) const
{
    return Kratos::make_shared<PotentialToCompressibleNavierStokesOperation>(rModel, Parameters);
}

const Parameters PotentialToCompressibleNavierStokesOperation::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"({
        "origin_model_part"       : "",
        "destination_model_part"  : "",
        "reference_temperature"   : 273,
        "compute_nodal_velocities": true
    })");
    return default_parameters;}

void PotentialToCompressibleNavierStokesOperation::Execute()
{
    KRATOS_TRY;

    const std::string origin_model_part_name = mParameters["origin_model_part"].GetString();
    const std::string destination_model_part_name = mParameters["destination_model_part"].GetString();
    const double& reference_temperature = mParameters["reference_temperature"].GetDouble();
    const double& compute_nodal_velocities = mParameters["compute_nodal_velocities"].GetBool();
    
    // Saving the modelparts
    auto& mOriginModelPart = mpModel->GetModelPart(origin_model_part_name);
    auto& mDestinationModelPart = mpModel->GetModelPart(destination_model_part_name);

    const int n_orig_nodes = mOriginModelPart.NumberOfNodes();
    const int n_dest_nodes = mDestinationModelPart.NumberOfNodes();

    // Check number of nodes
    KRATOS_ERROR_IF_NOT(n_orig_nodes == n_dest_nodes)
        << "Origin and destination model parts have different number of nodes."
        << "\n\t- Number of origin nodes: " << n_orig_nodes
        << "\n\t- Number of destination nodes: " << n_dest_nodes << std::endl;

    // Getting the required variables
    const double& gamma = mOriginModelPart.GetProcessInfo().GetValue(HEAT_CAPACITY_RATIO);
    const double& sound_velocity = mOriginModelPart.GetProcessInfo().GetValue(SOUND_VELOCITY);
    const double& free_stream_rho = mOriginModelPart.GetProcessInfo().GetValue(FREE_STREAM_DENSITY);
    const double& free_stream_mach = mOriginModelPart.GetProcessInfo().GetValue(FREE_STREAM_MACH);
    const double& specific_heat = (std::pow(sound_velocity,2) / (gamma * reference_temperature)) / (gamma - 1.0);
    
    // Calculate de nodal velocities
    if (compute_nodal_velocities){
        // Construct the ComputeNodalPotentialFlowVelocityProcess
        const std::vector<std::string> variable_array = {"VELOCITY"};
        ComputeNodalValueProcess ComputeNodalValueProcess(mOriginModelPart, variable_array);

        // Execute the ComputeNodalPotentialFlowVelocityProcess
        ComputeNodalValueProcess.Execute();}

    // Transfer the potential flow values as initial condition for the compressible problem
    IndexPartition<std::size_t>(n_orig_nodes).for_each([&](std::size_t index)
                                                       {
        auto it_dest_node = mDestinationModelPart.NodesBegin() + index;
        const auto it_orig_node = mOriginModelPart.NodesBegin() + index;

        // Getting potential flow velocities
        const auto &velocity = it_orig_node->GetValue(VELOCITY);
        
        // Calculate the conservative variables
        const double& velocity_norm_2 = velocity[0] * velocity[0] + velocity[1] * velocity[1];
        const double& velocity_norm = std::pow(velocity_norm_2,0.5);
        const double& mach = velocity_norm / sound_velocity;
        const double& num = 1.0 + 0.5 * (gamma - 1.0) * std::pow(free_stream_mach,2);
        const double& den = 1.0 + 0.5 * (gamma - 1.0) * std::pow(mach,2);
        const double& density = free_stream_rho * std::pow((num / den),(1.0 / (gamma - 1.0)));
        const double& energy = specific_heat * reference_temperature + 0.5 * velocity_norm_2;
        
        // Setting the Navier Stokes initial condition
        it_dest_node->FastGetSolutionStepValue(DENSITY, 0) = density; 
        it_dest_node->FastGetSolutionStepValue(MOMENTUM, 0) = density*velocity;
        it_dest_node->FastGetSolutionStepValue(TOTAL_ENERGY, 0) = density*energy; });

    KRATOS_CATCH("")
}

}
