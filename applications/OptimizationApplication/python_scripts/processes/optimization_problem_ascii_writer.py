from datetime import datetime
from typing import Union, Any

import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy import ExecutionPolicy
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem

def Factory(_: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ExecutionPolicy:
    if not parameters.Has("settings"):
        raise RuntimeError(f"IndependentAnalysisExecutionPolicy instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return OptimizationProblemAsciiWriter(parameters["settings"], optimization_problem)

class Header:
    def __init__(self, header_name: str, value: Any, format_info: dict, is_right_aligned: bool):
        header_name = header_name.strip()
        header_length = len(header_name)

        if isinstance(value, bool):
            value_length = 3
            value_format_post_fix = ""
            self.__value_converter = lambda x: "yes" if x else " no"
        elif isinstance(value, int):
            value_length = len(("{:" + str(format_info[type(value)]) + "d}").format(value))
            value_format_post_fix = "d"
            self.__value_converter = lambda x: int(x)
        elif isinstance(value, float):
            value_length = len(("{:0." + str(format_info[type(value)]) + "e}").format(value))
            value_format_post_fix = f".{format_info[type(value)]}e"
            self.__value_converter = lambda x: float(x)
        else:
            value_length = len(str(value))
            value_format_post_fix = "s"
            self.__value_converter = lambda x: str(x)

        if is_right_aligned:
            alignement_text = ">"
        else:
            alignement_text = "<"

        if header_length > value_length:
            self.__header_name = header_name
            self.__value_format = "{:" + alignement_text + str(header_length) + value_format_post_fix + "}"
        else:
            self.__header_name = ("{:" + alignement_text + str(value_length) + "s}").format(header_name)
            self.__value_format = "{:" + alignement_text + str(value_length) + value_format_post_fix + "}"

    def GetHeaderName(self) -> str:
        return self.__header_name

    def GetValueStr(self, value: Any) -> str:
        return self.__value_format.format(self.__value_converter(value))

class OptimizationProblemAsciiWriter(Kratos.OutputProcess):
    def GetDefaultParameters(self):
        return Kratos.Parameters(
            """
            {
                "output_file_name"         : "SPECIFY_OUTPUT_FILE_NAME",
                "write_kratos_version"     : true,
                "write_time_stamp"         : true,
                "write_initial_values"     : true,
                "list_of_output_components": ["all"],
                "format_info": {
                    "int_length"     : 7,
                    "float_precision": 9
                }
            }
            """
        )

    def __init__(self, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        Kratos.OutputProcess.__init__(self)

        self.optimization_problem = optimization_problem
        parameters.RecursivelyValidateAndAssignDefaults(self.GetDefaultParameters())

        self.output_file_name = parameters["output_file_name"].GetString()
        if not self.output_file_name.endswith(".csv"):
            self.output_file_name += ".csv"

        self.write_kratos_version = parameters["write_kratos_version"].GetBool()
        self.write_time_stamp = parameters["write_time_stamp"].GetBool()
        self.write_initial_values = parameters["write_initial_values"].GetBool()

        self.format_info = {
            int  : parameters["format_info"]["int_length"].GetInt(),
            float: parameters["format_info"]["float_precision"].GetInt()
        }

        self.list_of_components: 'list[Union[str, ResponseFunction, Control, ExecutionPolicy]]' = []
        list_of_component_names = parameters["list_of_output_components"].GetStringArray()
        if len(list_of_component_names) == 1 and list_of_component_names[0] == "all":
            # write all the data components
            # first add the responses
            for response_function in self.optimization_problem.GetListOfResponses():
                self.list_of_components.append(response_function)

            # then add controls
            for control in self.optimization_problem.GetListOfControls():
                self.list_of_components.append(control)

            # then add execution policies
            for execution_policy in self.optimization_problem.GetListOfExecutionPolicies():
                self.list_of_components.append(execution_policy)

            # now add the algorithm
            self.list_of_components.append("algorithm")
        else:
            for component_name in list_of_component_names:
                component_data = component_name.split(".")
                if len(component_data) == 2:
                    component_type = component_data[0]
                    component_name = component_data[1]
                    if component_type == "response_function":
                        self.list_of_components.append(self.optimization_problem.GetResponse(component_name))
                    elif component_type == "control":
                        self.list_of_components.append(self.optimization_problem.GetControl(component_name))
                    elif component_type == "execution_policy":
                        self.list_of_components.append(self.optimization_problem.GetExecutionPolicy(component_name))
                    else:
                        raise RuntimeError(f"Unsupported component type provided with component name string = \"{component_name}\". Supported component types are:\n\tresponse_function\n\tcontrol\n\texecution_policy\n\talgorithm")
                elif len(component_data) == 1:
                    self.list_of_components.append(component_name)
                else:
                    raise RuntimeError(f"Unsupported component type provided with component name string = \"{component_name}\". Supported component types are:\n\tresponse_function\n\tcontrol\n\texecution_policy\n\talgorithm")

        self.list_of_headers: 'list[tuple[Any, dict[str, Header]]]' = []
        self.initialized_headers = False

    def PrintOutput(self) -> None:
        if not self.initialized_headers:
            # now get the buffered data headers
            self.list_of_headers = self._GetHeaders(lambda x: x.GetBufferedData(), True)
            # write the ehader information
            self._WriteHeaders()
            self.initialized_headers = True

        if self._IsWritingProcess():
            # now write step data
            with open(self.output_file_name, "a") as file_output:
                # write the step
                file_output.write("{:>7d}".format(self.optimization_problem.GetStep()))

                # wrtie the values
                for component, header_info_dict in self.list_of_headers:
                    componend_data_view = ComponentDataView(component, self.optimization_problem)
                    buffered_dict = componend_data_view.GetBufferedData()
                    for k, header in header_info_dict.items():
                        file_output.write(", " + header.GetValueStr(buffered_dict[k]))

                file_output.write("\n")

    def ExecuteFinalize(self):
        if self._IsWritingProcess():
            with open(self.output_file_name, "a") as file_output:
                file_output.write("# End of file")

    def _IsWritingProcess(self):
        return True

    def _WriteHeaders(self):
        if (self._IsWritingProcess()):
            kratos_version = "not_given"
            if (self.write_kratos_version):
                kratos_version = str(Kratos.KratosGlobals.Kernel.Version())

            time_stamp = "not_specified"
            if (self.write_time_stamp):
                time_stamp = str(datetime.now())

            msg_header = ""
            msg_header += "# Optimization probelm ascii output\n"
            msg_header += f"# Kratos version: {kratos_version}\n"
            msg_header += f"# Timestamp     : {time_stamp}\n"
            msg_header += "# -----------------------------------------------\n"

            if self.write_initial_values:
                msg_header += "# --------------- Initial values ----------------\n"

                initial_headers = self._GetHeaders(lambda x: x.GetUnBufferedData(), False)
                # now write the initial value container data
                for component, header_info_dict in initial_headers:
                    componend_data_view = ComponentDataView(component, self.optimization_problem)
                    buffered_dict = componend_data_view.GetUnBufferedData()
                    if isinstance(component, str):
                        component_name = component
                    else:
                        component_name = component.GetName()
                    msg_header += "# \t" + component_name + ":\n"
                    for k, header in header_info_dict.items():
                        component_name_header = header.GetHeaderName()[len(component_name)+1:]
                        msg_header += "# \t\t" + component_name_header + ": " + header.GetValueStr(buffered_dict[k]) + "\n"

                msg_header += "# ------------ End of initial values ------------\n"
                msg_header += "# -----------------------------------------------\n"

            msg_header += "# ------------ Start of step values -------------\n"
            msg_header += "# Headers:\n"

            msg_header += "#  STEP"

            for _, header_info_dict in self.list_of_headers:
                for header in header_info_dict.values():
                    msg_header += ", " + header.GetHeaderName()

            msg_header += "\n"

            # write the header
            with open(self.output_file_name, "w") as file_output:
                file_output.write(msg_header)

    def _GetHeaders(self, dict_getter_method, is_right_aligned: bool) ->  'list[tuple[Any, dict[str, Header]]]':
        list_of_headers: 'list[tuple[Any, dict[str, Header]]]' = []
        for component in self.list_of_components:
            componend_data_view = ComponentDataView(component, self.optimization_problem)
            values_map = dict_getter_method(componend_data_view).GetMap()
            header_info_dict: 'dict[str, Header]' = {}
            if isinstance(component, str):
                component_name = component
            else:
                component_name = component.GetName()
            for k, v in values_map.items():
                if isinstance(v, (bool, int, float, str)):
                    header_name = component_name + ":" + k[k.rfind("/") + 1:]
                    if header_name in [header.GetHeaderName().strip() for header in header_info_dict.values()]:
                        Kratos.Logger.PrintWarning(self.__class__.__name__, "Second value with same header name = \"" + header_name + "\" found.")
                    header_info_dict[k] = Header(header_name, v, self.format_info, is_right_aligned)
            list_of_headers.append([component, header_info_dict])
        return list_of_headers