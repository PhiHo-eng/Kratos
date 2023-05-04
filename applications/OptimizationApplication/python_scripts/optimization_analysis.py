import KratosMultiphysics as Kratos
from KratosMultiphysics.process_factory import KratosProcessFactory
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm import Algorithm
from KratosMultiphysics.OptimizationApplication.model_part_controllers.model_part_controller import ModelPartController
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_decorator import ExecutionPolicyDecorator
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import OptimizationComponentFactory
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import CallOnAll
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem

class OptimizationAnalysis:
    @classmethod
    def GetDefaultParameters(cls) -> Kratos.Parameters:
        return Kratos.Parameters("""{
            "model_parts"       : [],
            "analyses"          : [],
            "responses"         : [],
            "controls"          : [],
            "algorithm_settings": {},
            "processes"  : {
                "kratos_processes"           : {},
                "optimization_data_processes": {}
            }
        }""")

    def __init__(self, model: Kratos.Model, project_parameters: Kratos.Parameters):
        self.model = model
        self.project_parameters = project_parameters
        self.project_parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.optimization_problem = OptimizationProblem(self.project_parameters["problem_data"]["echo_level"].GetInt())

        self.__list_of_model_part_controllers: 'list[ModelPartController]' = []
        self.__algorithm: Algorithm = None

        self._CreateModelPartControllers()
        self._CreateAnalyses()
        self._CreateControls()
        self._CreateResponses()
        self._CreateAlgorithm()
        self._CreateProcesses()

    def Initialize(self):
        CallOnAll(self.__list_of_model_part_controllers, ModelPartController.ImportModelPart)
        CallOnAll(self.__list_of_model_part_controllers, ModelPartController.Initialize)
        for process_type in self.__algorithm.GetProcessesOrder():
            CallOnAll(self.optimization_problem.GetListOfProcesses(process_type), Kratos.Process.ExecuteInitialize)
        CallOnAll(self.optimization_problem.GetListOfExecutionPolicies(), ExecutionPolicyDecorator.Initialize)
        CallOnAll(self.optimization_problem.GetListOfControls(), Control.Initialize)
        CallOnAll(self.optimization_problem.GetListOfResponses(), ResponseFunction.Initialize)

        self.__algorithm.Initialize()

    def Check(self):
        for process_type in self.__algorithm.GetProcessesOrder():
            CallOnAll(self.optimization_problem.GetListOfProcesses(process_type), Kratos.Process.Check)
        CallOnAll(self.optimization_problem.GetListOfExecutionPolicies(), ExecutionPolicyDecorator.Check)
        CallOnAll(self.optimization_problem.GetListOfControls(), Control.Check)
        CallOnAll(self.optimization_problem.GetListOfResponses(), ResponseFunction.Check)

        self.__algorithm.Check()

    def Finalize(self):
        self.__algorithm.Finalize()

        CallOnAll(self.__list_of_model_part_controllers, ModelPartController.Finalize)
        for process_type in self.__algorithm.GetProcessesOrder():
            CallOnAll(self.optimization_problem.GetListOfProcesses(process_type), Kratos.Process.ExecuteFinalize)
        CallOnAll(self.optimization_problem.GetListOfExecutionPolicies(), ExecutionPolicyDecorator.Finalize)
        CallOnAll(self.optimization_problem.GetListOfControls(), Control.Finalize)
        CallOnAll(self.optimization_problem.GetListOfResponses(), ResponseFunction.Finalize)

    def Run(self):
        self.Initialize()
        self.Check()
        self.__algorithm.SolveOptimizationProblem()
        self.Finalize()

    def _CreateModelPartControllers(self):
        default_model_part_controller_settings = Kratos.Parameters("""{
            "python_module": "mdpa_model_part_controller",
            "kratos_module": "KratosMultiphysics.OptimizationApplication.model_part_controllers",
            "Parameters"   : {}
        }""")
        # assign defaults to all model part controllers
        CallOnAll(self.project_parameters["model_parts"], Kratos.Parameters.AddMissingParameters, default_model_part_controller_settings)

        factory = KratosProcessFactory(self.model)
        self.__list_of_model_part_controllers: 'list[ModelPartController]' = factory.ConstructListOfProcesses(self.project_parameters["model_parts"])

    def _CreateAnalyses(self):
        for analyses_settings in self.project_parameters["analyses"]:
            execution_policy = ExecutionPolicyDecorator(self.model, analyses_settings, self.optimization_problem)
            self.optimization_problem.AddComponent(execution_policy)

    def _CreateResponses(self):
        default_settings = Kratos.Parameters("""{
            "name"         : "",
            "python_module": "",
            "kratos_module": "KratosMultiphysics.OptimizationApplication.responses",
            "Parameters"   : {}
        }""")
        for response_settings in self.project_parameters["responses"]:
            response_settings.ValidateAndAssignDefaults(default_settings)
            response_function: ResponseFunction = OptimizationComponentFactory(self.model, response_settings, self.optimization_problem)
            self.optimization_problem.AddComponent(response_function)

    def _CreateControls(self):
        default_settings = Kratos.Parameters("""{
            "name"          : "",
            "python_module" : "",
            "kratos_module" : "KratosMultiphysics.OptimizationApplication.controls",
            "Parameters"    : {}
        }""")
        for control_settings in self.project_parameters["controls"]:
            control_settings.ValidateAndAssignDefaults(default_settings)
            control = OptimizationComponentFactory(self.model, control_settings, self.optimization_problem)
            self.optimization_problem.AddComponent(control)

    def _CreateProcesses(self):
        process_settings = self.project_parameters["processes"]
        process_settings.ValidateAndAssignDefaults(self.GetDefaultParameters()["processes"])

        kratos_processes = process_settings["kratos_processes"]
        optimization_data_processes = process_settings["optimization_data_processes"]

        factory = KratosProcessFactory(self.model)

        process_lists_order = self.__algorithm.GetProcessesOrder()
        optimization_data_process_default_settings = Kratos.Parameters("""{
            "python_module" : "",
            "kratos_module" : "KratosMultiphysics.OptimizationApplication.optimization_data_processes",
            "Parameters"    : {}
        }""")

        for process_type in process_lists_order:
            if kratos_processes.Has(process_type):
                for process in factory.ConstructListOfProcesses(kratos_processes[process_type]):
                    self.optimization_problem.AddProcess(process_type, process)
            if optimization_data_processes.Has(process_type):
                for process_settings in optimization_data_processes[process_type]:
                    process_settings.ValidateAndAssignDefaults(optimization_data_process_default_settings)
                    process = OptimizationComponentFactory(self.model, process_settings, self.optimization_problem)
                    self.optimization_problem.AddProcess(process_type, process)

    def _CreateAlgorithm(self):
        default_settings = Kratos.Parameters("""{
            "python_module" : "",
            "kratos_module" : "KratosMultiphysics.OptimizationApplication.algorithms",
            "Parameters"    : {}
        }""")
        algorithm_settings = self.project_parameters["algorithm_settings"]
        algorithm_settings.ValidateAndAssignDefaults(default_settings)
        self.__algorithm = OptimizationComponentFactory(self.model, algorithm_settings, self.optimization_problem)
