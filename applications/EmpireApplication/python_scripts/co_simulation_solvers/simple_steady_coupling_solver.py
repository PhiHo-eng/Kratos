# Importing the base class
from co_simulation_solvers.co_simulation_base_coupling_solver import CoSimulationBaseCouplingSolver
import KratosMultiphysics
# Other imports
from co_simulation_convergence_accelerators.co_simulation_convergence_accelerator_factory import CreateConvergenceAccelerator
from co_simulation_convergence_criteria.co_simulation_convergence_criteria_factory import CreateConvergenceCriteria
from co_simulation_tools import couplingsolverprint, red, green, cyan, bold


def CreateSolver(cosim_solver_settings, level):
    return SimpleSteadyCouplingSolver(cosim_solver_settings, level)

class SimpleSteadyCouplingSolver(CoSimulationBaseCouplingSolver):
    def __init__(self, cosim_solver_settings, level):
        if not len(cosim_solver_settings["solvers"]) == 2:
            raise Exception("Exactly two solvers have to be specified for the " + self.__class__.__name__ + "!")

        super(SimpleSteadyCouplingSolver, self).__init__(cosim_solver_settings, level)
        self.time = self.cosim_solver_settings["start_coupling_time"]

        self.convergence_accelerator = CreateConvergenceAccelerator(
            self.cosim_solver_settings["convergence_accelerator_settings"],
            self.solvers, self.lvl)
        self.convergence_accelerator.SetEchoLevel(self.echo_level)

        self.convergence_criteria = CreateConvergenceCriteria(
            self.cosim_solver_settings["convergence_criteria_settings"],
            self.solvers, self.lvl)
        self.convergence_criteria.SetEchoLevel(self.echo_level)

        self.num_coupling_iterations = self.cosim_solver_settings["num_coupling_iterations"]


    def Initialize(self):
        super(SimpleSteadyCouplingSolver, self).Initialize()
        self.convergence_accelerator.Initialize()
        self.convergence_criteria.Initialize()

    def Finalize(self):
        super(SimpleSteadyCouplingSolver, self).Finalize()
        self.convergence_accelerator.Finalize()
        self.convergence_criteria.Finalize()

    def InitializeSolutionStep(self):
        super(SimpleSteadyCouplingSolver, self).InitializeSolutionStep()
        self.convergence_accelerator.InitializeSolutionStep()
        self.convergence_criteria.InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        super(SimpleSteadyCouplingSolver, self).FinalizeSolutionStep()
        self.convergence_accelerator.FinalizeSolutionStep()
        self.convergence_criteria.FinalizeSolutionStep()

    def SolveSolutionStep(self):

        for k in range(self.num_coupling_iterations):
            if self.echo_level > 0:
                couplingsolverprint(self.lvl, self._Name(),
                                    cyan("Coupling iteration:"), bold(str(k+1)+" / " + str(self.num_coupling_iterations)))

            self.convergence_accelerator.InitializeNonLinearIteration()
            self.convergence_criteria.InitializeNonLinearIteration()

            for solver_name in self.solver_names:
                print('solver_name', solver_name)
                solver = self.solvers[solver_name]
                self._SynchronizeInputData(solver, solver_name)
                solver.SolveSolutionStep()
                self._SynchronizeOutputData(solver, solver_name)

            self.convergence_accelerator.FinalizeNonLinearIteration()
            self.convergence_criteria.FinalizeNonLinearIteration()

            if self.convergence_criteria.IsConverged():
                if self.echo_level > 0:
                    couplingsolverprint(self.lvl, self._Name(), green("### CONVERGENCE WAS ACHIEVED ###"))
                break
            else:
                self.convergence_accelerator.ComputeUpdate()

            if k+1 >= self.num_coupling_iterations and self.echo_level > 0:
                couplingsolverprint(self.lvl, self._Name(), red("XXX CONVERGENCE WAS NOT ACHIEVED XXX"))

    def _SynchronizeInputData(self, solver, solver_name):

        input_data_list = self.cosim_solver_details[solver_name]["input_data_list"]

        for input_data in input_data_list:
            from_solver = self.solvers[input_data["from_solver"]]
            data_name = input_data["data_name"]
            data_definition = from_solver.GetDataDefinition(data_name)
            data_settings = { "data_format" : data_definition["data_format"],
                            "data_name"   : data_name,
                            "io_settings" : input_data["io_settings"] }
            solver.ImportData(data_settings, from_solver)

    def _SynchronizeOutputData(self, solver, solver_name):
        output_data_list = self.cosim_solver_details[solver_name]["output_data_list"]

        for output_data in output_data_list:
            to_solver = self.solvers[output_data["to_solver"]]
            data_name = output_data["data_name"]
            data_definition = to_solver.GetDataDefinition(data_name)
            data_settings = { "data_format" : data_definition["data_format"],
                            "data_name"   : data_name,
                            "io_settings" : output_data["io_settings"] }
            solver.ExportData(data_settings, to_solver)
