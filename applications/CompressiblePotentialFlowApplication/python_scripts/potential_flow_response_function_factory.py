
from KratosMultiphysics.CompressiblePotentialFlowApplication import potential_flow_response

try:
    import KratosMultiphysics.MultilevelMonteCarloApplication
    import KratosMultiphysics.MappingApplication
    import xmc
    import exaqute
    import numpy as np
    from KratosMultiphysics.CompressiblePotentialFlowApplication import stochastic_potential_flow_response
    is_xmc_available = True
except Exception as err:
    print(err)
    stop
    is_xmc_available = False

def CreateResponseFunction(response_id, response_settings, model):
    response_type = response_settings["response_type"].GetString()

    if response_type == "adjoint_lift_potential_jump":
        return potential_flow_response.AdjointResponseFunction(response_id, response_settings, model)
    elif response_type == "embedded_adjoint_lift_potential_jump":
        return potential_flow_response.EmbeddedAdjointResponseFunction(response_id, response_settings, model)
    elif response_type == "angle_of_attack":
        return potential_flow_response.AngleOfAttackResponseFunction(response_id, response_settings, model)
    elif response_type == "chord_length":
        return potential_flow_response.ChordLengthResponseFunction(response_id, response_settings, model)
    elif response_type == "perimeter":
        return potential_flow_response.PerimeterResponseFunction(response_id, response_settings, model)
    elif response_type == "adjoint_lift_potential_jump_with_perimeter":
        return potential_flow_response.LiftPerimeterResponseFunction(response_id, response_settings, model)
    elif response_type == "stochastic_adjoint_lift_potential_jump":
        if is_xmc_available:
            return stochastic_potential_flow_response.AdjointResponseFunction(response_id, response_settings, model)
        else:
            err_msg = "XMC and its dependencies could not be imported. Please check applications/MultilevelMonteCarloApplication/README.md"
            err_msg += " for installation details. Please make sure that the CompressiblePotentialFlowApplication, MultilevelMonteCarloApplication"
            err_msg += " and the MappingApplication are compiled."
            raise ImportError(err_msg)
    else:
        raise NameError("The type of the following response function is not specified: "+ response_id +
                        ".\nAvailable types are: 'adjoint_lift_potential_jump', 'angle_of_attack', 'chord_length'," +
                        ",'perimeter' and 'adjoint_lift_potential_jump_with_perimeter'." )
