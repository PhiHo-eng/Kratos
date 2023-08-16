try:
    import nlopt
except ImportError:
    raise Exception("NLOPT python library is not available")

import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.algorithms.standardized_NLOPT_objective import StandardizedNLOPTObjective
from KratosMultiphysics.OptimizationApplication.algorithms.standardized_NLOPT_constraint import StandardizedNLOPTConstraint
from KratosMultiphysics.OptimizationApplication.controls.master_control import MasterControl
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm import Algorithm
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import DictLogger
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import CallOnAll


def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
    return NLOPTAlgorithms(model, parameters, optimization_problem)

class NLOPTAlgorithms(Algorithm):
    """
        A classical steepest descent algorithm to solve unconstrainted optimization problems.
    """

    @classmethod
    def GetDefaultParameters(cls):
        return Kratos.Parameters("""{
            "module"            : "KratosMultiphysics.OptimizationApplication.algorithms",
            "type"              : "PLEASE_PROVIDE_AN_ALGORITHM_CLASS_NAME",
            "objective"         : {},
            "constraints"         : [],
            "controls"          : [],
            "echo_level"        : 0,
            "NLOPT_settings"          : {
                "algorithm_name"      : "MMA",
                "verbosity"           : 0,
                "controls_lower_bound": "",
                "controls_upper_bound": "",
                "stopping_criteria"   : {
                   "objective_rel_tol": 1e-6,
                   "objective_abs_tol": 1e-6,
                   "controls_rel_tol": 1e-6,
                   "controls_abs_tol": 1e-6,
                   "maximum_function_evalualtion": 10
                },
                "algorithm_specific_settings"   : {}
            }
        }""")

    def __init__(self, model:Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        self.model = model
        self.parameters = parameters
        self._optimization_problem = optimization_problem

        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())

        # controls
        self.master_control = MasterControl()
        for control_name in parameters["controls"].GetStringArray():
            control = optimization_problem.GetControl(control_name)
            self.master_control.AddControl(control)
        self.__control_field = None

        # objective & constraints
        self.__objective = StandardizedNLOPTObjective(parameters["objective"], self.master_control, self._optimization_problem)
        self.__constraints = []
        for constraint_settings in parameters["constraints"]:
            self.__constraints.append(StandardizedNLOPTConstraint(constraint_settings, self.master_control, self._optimization_problem))

        # nlopt settings
        NLOPT_settings = parameters["NLOPT_settings"]
        NLOPT_settings.ValidateAndAssignDefaults(self.GetDefaultParameters()["NLOPT_settings"])

        # nlopt verbosity
        self.nlopt_verbosity = NLOPT_settings["verbosity"].GetInt()

        # upper and lower bounds
        self.nlopt_controls_lower_bound = NLOPT_settings["controls_lower_bound"].GetString()
        self.nlopt_controls_upper_bound = NLOPT_settings["controls_upper_bound"].GetString()

        # nlopt algorithm
        self.algorithm_name = NLOPT_settings["algorithm_name"].GetString()
        if(len(self.__constraints)==0):
            self.CheckOptimizerSupport(self.algorithm_name,None)
        else:
            for constraint in self.__constraints:
                if constraint.IsEqualityType():
                    self.CheckOptimizerSupport(self.algorithm_name,"equality")
                else:
                    self.CheckOptimizerSupport(self.algorithm_name,"inequality")

        # stopping
        self.stopping_criteria = NLOPT_settings["stopping_criteria"]
        self.stopping_criteria.ValidateAndAssignDefaults(self.GetDefaultParameters()["NLOPT_settings"]["stopping_criteria"])

        # alg specific settings
        self.opt_algorithm_specific_settings = NLOPT_settings["algorithm_specific_settings"]

    def GetMinimumBufferSize(self) -> int:
        return 2

    def Check(self):
        pass

    def Initialize(self):
        self.converged = False
        self.__objective.Initialize()
        self.__objective.Check()
        for constraint in self.__constraints:
            constraint.Initialize()
            constraint.Check()
        self.master_control.Initialize()
        self.__control_field = self.master_control.GetControlField()
        self.algorithm_data = ComponentDataView("algorithm", self._optimization_problem)

        # create nlopt optimizer
        self.x0 = self.__control_field.Evaluate()
        self.x0 = self.x0.reshape(-1)
        self.nlopt_optimizer = nlopt.opt(self.GetOptimizer(self.algorithm_name), self.x0.size)

        # set nlopt verbosity
        self.nlopt_optimizer.set_param("verbosity",self.nlopt_verbosity)

        # assign objectives and constarints
        self.nlopt_optimizer.set_min_objective(self.__objective.CalculateStandardizedValueAndGradients)
        for constraint in self.__constraints:
            if constraint.IsEqualityType():
                self.nlopt_optimizer.add_equality_constraint(lambda x,grad: constraint.CalculateStandardizedValueAndGradients(x,grad),1e-8)
            else:
                self.nlopt_optimizer.add_inequality_constraint(lambda x,grad: constraint.CalculateStandardizedValueAndGradients(x,grad),1e-8)

        # now set stopping criteria
        self.nlopt_optimizer.set_ftol_rel(self.stopping_criteria["objective_rel_tol"].GetDouble())
        self.nlopt_optimizer.set_ftol_abs(self.stopping_criteria["objective_abs_tol"].GetDouble())
        self.nlopt_optimizer.set_xtol_rel(self.stopping_criteria["controls_rel_tol"].GetDouble())
        self.nlopt_optimizer.set_xtol_abs(self.stopping_criteria["controls_abs_tol"].GetDouble())
        self.nlopt_optimizer.set_maxeval(self.stopping_criteria["maximum_function_evalualtion"].GetInt())

        # set bounds if exist
        if self.nlopt_controls_lower_bound != "":
            self.nlopt_optimizer.set_lower_bounds(float(self.nlopt_controls_lower_bound))
        if self.nlopt_controls_upper_bound != "":
            self.nlopt_optimizer.set_upper_bounds(float(self.nlopt_controls_upper_bound))

        # now add algorithm specific settings
        if self.opt_algorithm_specific_settings.Has("inner_maxeval"):
            self.nlopt_optimizer.set_param("inner_maxeval",self.opt_algorithm_specific_settings["inner_maxeval"].GetInt())
        if self.opt_algorithm_specific_settings.Has("dual_ftol_rel"):
            self.nlopt_optimizer.set_param("dual_ftol_rel",self.opt_algorithm_specific_settings["dual_ftol_rel"].GetDouble())
        if self.opt_algorithm_specific_settings.Has("initial_step"):
            self.nlopt_optimizer.set_initial_step(self.opt_algorithm_specific_settings["initial_step"].GetDouble())
        if self.opt_algorithm_specific_settings.Has("population"):
            self.nlopt_optimizer.set_population(self.opt_algorithm_specific_settings["population"].GetInt())
        if self.opt_algorithm_specific_settings.Has("number_of_stored_vectors"):
            self.nlopt_optimizer.set_vector_storage(self.opt_algorithm_specific_settings["number_of_stored_vectors"].GetInt())

    def Finalize(self):
        self.__objective.Finalize()
        for constraint in self.__constraints:
            constraint.Finalize()
        self.master_control.Finalize()

    def Solve(self):
        self.nlopt_optimizer.optimize(self.x0)
        CallOnAll(self._optimization_problem.GetListOfProcesses("output_processes"), Kratos.OutputProcess.PrintOutput)
        self.LogTermination()

    def LogTermination(self):
        stopval = self.nlopt_optimizer.get_stopval()
        maxtime = self.nlopt_optimizer.get_maxtime()
        nlopt_result_map = {
            nlopt.SUCCESS: 'Generic success return value.',
            nlopt.STOPVAL_REACHED: f'Optimization stopped because stopval ({stopval}) was reached.',
            nlopt.FTOL_REACHED: f'Optimization stopped because ftol_rel ({self.stopping_criteria["objective_rel_tol"].GetDouble()}) or ftol_abs ({self.stopping_criteria["objective_abs_tol"].GetDouble()}) was reached.',
            nlopt.XTOL_REACHED: f'Optimization stopped because xtol_rel ({self.stopping_criteria["controls_rel_tol"].GetDouble()}) or xtol_abs ({self.stopping_criteria["controls_abs_tol"].GetDouble()}) was reached.',
            nlopt.MAXEVAL_REACHED: f'Optimization stopped because maxeval ({self.stopping_criteria["maximum_function_evalualtion"].GetInt()}) was reached.',
            nlopt.MAXTIME_REACHED: f'Optimization stopped because maxtime ({maxtime}) was reached.',
            nlopt.FAILURE: 'Generic failure return value.',
            nlopt.INVALID_ARGS: 'Invalid arguments (e.g. lower bounds are bigger than upper bounds, an unknown algorithm was specified, etcetera).',
            nlopt.OUT_OF_MEMORY: 'Ran out of memory.',
            nlopt.ROUNDOFF_LIMITED: 'Halted because roundoff errors limited progress. (In this case, the optimization still typically returns a useful result.)',
            nlopt.FORCED_STOP: 'Halted because of a forced termination: the user called nlopt_force_stop(opt) on the optimization’s nlopt_opt object opt from the user’s objective function or constraints.'
        }
        info = {"Termination":nlopt_result_map[self.nlopt_optimizer.last_optimize_result()]}
        DictLogger(f"NLOPT-{self.algorithm_name}",info)

    def CheckOptimizerSupport(self, algorithm_name, constraint_type):
        # Unconstrained gradient-based optimizers
        unconstrained_optimizers = {"LBFGS","MMA","CCSAQ","TNEWTON_PRECOND_RESTART","TNEWTON_PRECOND","TNEWTON_RESTART"}

        # Constrained gradient-based optimizers with equality constraints
        equality_constrained_optimizers = {"SLSQP"}

        # Constrained gradient-based optimizers with inequality constraints
        inequality_constrained_optimizers = {"MMA","SLSQP","CCSAQ"}

        if constraint_type == "equality":
            if algorithm_name not in equality_constrained_optimizers:
                raise ValueError(f"The algorithm {algorithm_name} does not support equality constraints, support {equality_constrained_optimizers}.")
        elif constraint_type == "inequality":
            if algorithm_name not in inequality_constrained_optimizers:
                raise ValueError(f"The algorithm {algorithm_name} does not support inequality constraints, support {inequality_constrained_optimizers}.")
        else:
            if algorithm_name not in unconstrained_optimizers:
                raise ValueError(f"The algorithm {algorithm_name} does not support unconstrained optimization, support {unconstrained_optimizers}.")
        return True

    def GetOptimizer(self,name):
        nlopt_algorithm_mapping = {
            "MMA": nlopt.LD_MMA,
            "SLSQP": nlopt.LD_SLSQP,
            "LBFGS": nlopt.LD_LBFGS,
            "CCSAQ": nlopt.LD_CCSAQ,
            "TNEWTON_PRECOND_RESTART": nlopt.LD_TNEWTON_PRECOND_RESTART,
            "TNEWTON_PRECOND": nlopt.LD_TNEWTON_PRECOND,
            "TNEWTON_RESTART": nlopt.LD_TNEWTON_RESTART
        }
        return nlopt_algorithm_mapping[name]
