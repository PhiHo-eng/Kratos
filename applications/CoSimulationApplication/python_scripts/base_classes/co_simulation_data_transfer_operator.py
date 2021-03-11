# Importing the Kratos Library
import KratosMultiphysics as KM

class CoSimulationDataTransferOperator:
    """Baseclass for the data transfer operators used for CoSimulation
    It transfers data from one interface to another. This can e.g. be mapping or a copy of values.
    """
    def __init__(self, settings):
        self.settings = settings
        default_parameters = self._GetDefaultParameters()
        self.settings.ValidateAndAssignDefaults(default_parameters)
        self.settings["update_interface"].ValidateAndAssignDefaults(default_parameters["update_interface"])

        self.echo_level = self.settings["echo_level"].GetInt()
        self.__checked_combinations = []
        self.update_interface_execution_points = self.settings["update_interface"]["execution_points"].GetStringArray()

        available_execution_points = [
            "initialize_solution_step",
            "finalize_solution_step",
            "initialize_non_linear_iteration",
            "finalize_non_linear_iteration"
        ]

        for execution_point in self.update_interface_execution_points:
            if execution_point not in available_execution_points:
                err_msg  = 'Execution point "{}" is not available, only the following options are available:\n    '.format(execution_point)
                err_msg += "\n    ".join(available_execution_points)
                raise Exception(err_msg)

    def InitializeSolutionStep(self):
        if "initialize_solution_step" in self.update_interface_execution_points:
            self._UpdateInterface(self.settings["update_interface"]["settings"])

    def FinalizeSolutionStep(self):
        if "finalize_solution_step" in self.update_interface_execution_points:
            self._UpdateInterface(self.settings["update_interface"]["settings"])

    def InitializeNonLinearIteration(self):
        if "initialize_non_linear_iteration" in self.update_interface_execution_points:
            self._UpdateInterface(self.settings["update_interface"]["settings"])

    def FinalizeNonLinearIteration(self):
        if "finalize_non_linear_iteration" in self.update_interface_execution_points:
            self._UpdateInterface(self.settings["update_interface"]["settings"])

    def TransferData(self, from_solver_data, to_solver_data, transfer_options):
        # 1. Check if specified transfer options are available
        self._CheckAvailabilityTransferOptions(transfer_options)

        # 2. Perform check (only if it has not been done before in this combination)
        identifier_from_solver_data = from_solver_data.solver_name + "." + from_solver_data.model_part_name
        identifier_to_solver_data   = to_solver_data.solver_name   + "." + to_solver_data.model_part_name

        identifier_tuple = (identifier_from_solver_data, identifier_to_solver_data)
        if not identifier_tuple in self.__checked_combinations:
            self.__checked_combinations.append(identifier_tuple)
            self._Check(from_solver_data, to_solver_data)

        # 3. Perform data transfer
        self._ExecuteTransferData(from_solver_data, to_solver_data, transfer_options)

    def _ExecuteTransferData(self, from_solver_data, to_solver_data, transfer_options):
        raise NotImplementedError("This function has to be implemented in the derived class!")

    def _UpdateInterface(self, settings):
        raise NotImplementedError("This function is not implemented for {}".format(self._ClassName()))

    def _Check(self, from_solver_data, to_solver_data):
        # this can be implemented in derived classes if necessary
        # the purpose is to check only once each combination of data
        # this mechanism is necessary since the data-transfer operators can be used for different
        # combinations of data on the fly, i.e. they cannot be checked after the Initialization
        pass

    @classmethod
    def _ClassName(cls):
        return cls.__name__

    @classmethod
    def _GetListAvailableTransferOptions(cls):
        raise NotImplementedError("This function has to be implemented in the derived class!")

    @classmethod
    def _CheckAvailabilityTransferOptions(cls, transfer_options):
        avail_options = cls._GetListAvailableTransferOptions()
        for option_name in transfer_options.GetStringArray():
            if not option_name in avail_options:
                err_msg  = 'transfer option "{}" not recognized for "{}"!\n'.format(option_name, cls._ClassName())
                err_msg += 'Available options: "{}"'.format('", "'.join(avail_options))
                raise Exception(err_msg)

    @classmethod
    def _GetDefaultParameters(cls):
        return KM.Parameters("""{
            "type"       : "UNSPECIFIED",
            "echo_level" : 0,
            "update_interface" : {
                "execution_points" : [],
                "settings" : {}
            }
        }""")
