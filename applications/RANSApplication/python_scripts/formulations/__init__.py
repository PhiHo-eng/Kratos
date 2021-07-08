__all__ = ["Factory"]

import KratosMultiphysics as Kratos

# stabilization validation formulations
from KratosMultiphysics.RANSApplication.formulations.stabilization_validation.body_force_governed_cdr_rans_formulation import BodyForceGovernedCDRRansFormulation
from KratosMultiphysics.RANSApplication.formulations.stabilization_validation.circular_convection_rans_formulation import CircularConvectionRansFormulation

# flow solver formulations
from KratosMultiphysics.RANSApplication.formulations.incompressible_potential_flow import IncompressiblePotentialFlowRansFormulation
from KratosMultiphysics.RANSApplication.formulations.monolithic_vms.monolithic_velocity_pressure_rans_formulation import MonolithicVelocityPressureRansFormulation
from KratosMultiphysics.RANSApplication.formulations.fractional_step.fractional_step_velocity_pressure_rans_formulation import FractionalStepVelocityPressureRansFormulation

# turbulence model formulations
### k-epsilon formulations
from KratosMultiphysics.RANSApplication.formulations.turbulence_models.k_epsilon_rans_formulation import KEpsilonRansFormulation
from KratosMultiphysics.RANSApplication.formulations.monolithic_vms.monolithic_k_epsilon_rans_formulation import MonolithicKEpsilonRansFormulation
from KratosMultiphysics.RANSApplication.formulations.fractional_step.fractional_step_k_epsilon_rans_formulation import FractionalStepKEpsilonRansFormulation

### k-omega formulations
from KratosMultiphysics.RANSApplication.formulations.turbulence_models.k_omega_rans_formulation import KOmegaRansFormulation
from KratosMultiphysics.RANSApplication.formulations.monolithic_vms.monolithic_k_omega_rans_formulation import MonolithicKOmegaRansFormulation
from KratosMultiphysics.RANSApplication.formulations.fractional_step.fractional_step_k_omega_rans_formulation import FractionalStepKOmegaRansFormulation

### k-omega-sst formulations
from KratosMultiphysics.RANSApplication.formulations.turbulence_models.k_omega_sst_rans_formulation import KOmegaSSTRansFormulation
from KratosMultiphysics.RANSApplication.formulations.monolithic_vms.monolithic_k_omega_sst_rans_formulation import MonolithicKOmegaSSTRansFormulation
from KratosMultiphysics.RANSApplication.formulations.fractional_step.fractional_step_k_omega_sst_rans_formulation import FractionalStepKOmegaSSTRansFormulation

def Factory(model_part, settings):
    formulation_name = settings["formulation_name"].GetString()
    formulations_list = [
        # adding flow solvers including turbulence solvers
        ["incompressible_potential_flow", IncompressiblePotentialFlowRansFormulation],
        ["monolithic", MonolithicVelocityPressureRansFormulation],
        ["fractional_step", FractionalStepVelocityPressureRansFormulation],
        ["monolithic_k_epsilon", MonolithicKEpsilonRansFormulation],
        ["fractional_step_k_epsilon", FractionalStepKEpsilonRansFormulation],
        ["monolithic_k_omega", MonolithicKOmegaRansFormulation],
        ["fractional_step_k_omega", FractionalStepKOmegaRansFormulation],
        ["monolithic_k_omega_sst", MonolithicKOmegaSSTRansFormulation],
        ["fractional_step_k_omega_sst", FractionalStepKOmegaSSTRansFormulation],

        # adding only the turbulence models, flow should be given as a solved quantity for following formulations
        ["k_epsilon", KEpsilonRansFormulation],
        ["k_omega", KOmegaRansFormulation],
        ["k_omega_sst", KOmegaSSTRansFormulation],

        # adding stabilization validation formulations
        ["circular_convection", CircularConvectionRansFormulation],
        ["body_force_governed_cdr", BodyForceGovernedCDRRansFormulation]
    ]

    formulation_names_list = [
        formulations_list[i][0] for i in range(len(formulations_list))
    ]
    formulation_list = [
        formulations_list[i][1] for i in range(len(formulations_list))
    ]

    if (formulation_name not in formulation_names_list):
        msg = "Unknown formulation_name = \"" + formulation_name + "\". \nFollowing formulations are allowed:\n    "
        msg += "\n    ".join(sorted(formulation_names_list))
        raise Exception(msg + "\n")

    current_formulation = formulation_list[formulation_names_list.index(
        formulation_name)](model_part, settings)

    Kratos.Logger.PrintInfo("RansFormulationFactory",
                            "Created " + formulation_name + " formulation.")

    return current_formulation
