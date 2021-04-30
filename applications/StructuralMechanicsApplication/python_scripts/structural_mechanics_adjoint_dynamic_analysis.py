from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics
from KratosMultiphysics.StructuralMechanicsApplication import python_solvers_wrapper_structural as structural_solvers

# Importing the base class
from KratosMultiphysics.analysis_stage import AnalysisStage

class StructuralMechanicsAdjointDynamicAnalysis(AnalysisStage):
    """
    This class is the main-script of the StructuralMechanicsApplication put in a class

    It can be imported and used as "black-box"
    """
    def __init__(self, model, project_parameters):
        # Making sure that older cases still work by properly initalizing the parameters
        solver_settings = project_parameters["solver_settings"]

        if solver_settings.Has("domain_size") and project_parameters["problem_data"].Has("domain_size"):
            raise Exception("StructuralMechanicsAnalysis: " + '"domain_size" defined both in "problem_data" and "solver_settings"!')

        if solver_settings.Has("model_part_name") and project_parameters["problem_data"].Has("model_part_name"):
            raise Exception("StructuralMechanicsAnalysis: " + '"model_part_name" defined both in problem_data" and "solver_settings"!')

        if solver_settings.Has("time_stepping") and project_parameters["problem_data"].Has("time_Step"):
            raise Exception("StructuralMechanicsAnalysis: " + '"time_stepping" defined both in "problem_data" and "solver_settings"!')

        if not solver_settings.Has("time_stepping"):
            raise Exception("StructuralMechanicsAnalysis: Using the old way to pass the time_step, this was removed!")

        if not solver_settings.Has("domain_size"):
            raise Exception("StructuralMechanicsAnalysis: Using the old way to pass the domain_size, this was removed!")

        if not project_parameters["problem_data"].Has("end_time"):
            project_parameters["problem_data"].AddEmptyValue("end_time")
            project_parameters["problem_data"]["end_time"].SetDouble( \
                            project_parameters["problem_data"]["start_step"].GetDouble() + \
                            project_parameters["problem_data"]["nsteps"].GetInt()*solver_settings["time_stepping"]["time_step"].GetDouble()
                        )
        
        if not project_parameters["problem_data"].Has("start_time"):
            project_parameters["problem_data"].AddEmptyValue("start_time")
            project_parameters["problem_data"]["start_time"].SetDouble( \
                            project_parameters["problem_data"]["start_step"].GetDouble() \
                            )

        self.number_of_steps = project_parameters["problem_data"]["nsteps"].GetInt()

        if not solver_settings["response_function_settings"].Has("time_domain"):
            solver_settings["response_function_settings"].AddEmptyValue("time_domain")
            solver_settings["response_function_settings"]["time_domain"].SetDouble( \
                            project_parameters["problem_data"]["start_time"].GetDouble() - \
                            project_parameters["problem_data"]["end_time"].GetDouble()
                            )
        # TODO move this
    
        # Detect is a contact problem
        # NOTE: We have a special treatment for contact problems due to the way the convergence info is printed (in a table). Not doing this will provoque that the table is discontinous (and not fancy and eye-candy)
        solver_settings = project_parameters["solver_settings"]
        self.contact_problem = solver_settings.Has("contact_settings") or solver_settings.Has("mpc_contact_settings")

        super().__init__(model, project_parameters)

    def Initialize(self):
        """ Initializing the Analysis """

        # checking if "USE_CONSISTENT_MASS_MATRIX" is used in the materials file and correct it
        mat_import_settings = self.project_parameters["solver_settings"]["material_import_settings"]
        if mat_import_settings.Has("materials_filename"):
            materials_filename = mat_import_settings["materials_filename"].GetString()
            if materials_filename != "": # user specified the materials
                # check if old variables are used
                with open(materials_filename,'r') as mat_file:
                    materials = KratosMultiphysics.Parameters(mat_file.read())
                    for mat in materials["properties"]:
                        if mat.Has("Material"):
                            mat_spec = mat["Material"]
                            if mat_spec.Has("Variables"):
                                variables = mat_spec["Variables"]
                                for var_name in variables.keys():
                                    if var_name == "USE_CONSISTENT_MASS_MATRIX":
                                        raise Exception('Variable "USE_CONSISTENT_MASS_MATRIX" found. This was replaced by "COMPUTE_LUMPED_MASS_MATRIX"! Please adapt your input')
                                        
        super().Initialize()

        # dummy time step to correctly calculate DELTA_TIME
        self._GetSolver().main_model_part.CloneTimeStep(self.time)

        # In case of contact problem
        if self.contact_problem:
            self._GetSolver().SetEchoLevel(self.echo_level)
            # To avoid many prints
            if self.echo_level == 0:
                KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

    def RunSolutionLoop(self):
        """Note that the adjoint problem is solved in reverse time
        """
        for step in range(self.number_of_steps):
            self.time = self._GetSolver().AdvanceInTime(self.time)
            self.InitializeSolutionStep()
            self._GetSolver().Predict()
            self._GetSolver().SolveSolutionStep()
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()

    def OutputSolutionStep(self):
        """This function printed / writes output files after the solution of a step
        """

        # In case of contact problem
        if self.contact_problem:
            # First we check if one of the output processes will print output in this step this is done to save computation in case none of them will print
            is_output_step = False
            for output_process in self._GetListOfOutputProcesses():
                if output_process.IsOutputStep():
                    is_output_step = True
                    break

            if is_output_step:
                # Informing the output will be created
                KratosMultiphysics.Logger.PrintWarning("StructuralMechanicsAnalysis", "STEP: ", self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.STEP])
                KratosMultiphysics.Logger.PrintWarning("StructuralMechanicsAnalysis", "TIME: ", self.time)

        # Creating output
        super().OutputSolutionStep()

    #### Internal functions ####
    def _CreateSolver(self):
        """ Create the Solver (and create and import the ModelPart if it is not alread in the model) """
        ## Solver construction
        return structural_solvers.CreateSolver(self.model, self.project_parameters)

    def _CreateProcesses(self, parameter_name, initialization_order):
        """Create a list of Processes
        This method is TEMPORARY to not break existing code
        It will be removed in the future
        """
        list_of_processes = super()._CreateProcesses(parameter_name, initialization_order)

        if parameter_name == "processes":
            processes_block_names = ["constraints_process_list", "loads_process_list", "list_other_processes", "json_output_process",
                "json_check_process", "check_analytic_results_process", "contact_process_list"]
            if len(list_of_processes) == 0: # Processes are given in the old format (or no processes are specified)
                for process_name in processes_block_names:
                    if self.project_parameters.Has(process_name):
                        info_msg  = "Using the old way to create the processes, this was removed!\n"
                        info_msg += "Refer to \"https://github.com/KratosMultiphysics/Kratos/wiki/Common-"
                        info_msg += "Python-Interface-of-Applications-for-Users#analysisstage-usage\" "
                        info_msg += "for a description of the new format"
                        raise Exception("StructuralMechanicsAnalysis: " + info_msg)

            else: # Processes are given in the new format
                for process_name in processes_block_names:
                    if self.project_parameters.Has(process_name):
                        raise Exception("Mixing of process initialization is not allowed!")
        elif parameter_name == "output_processes":
            if self.project_parameters.Has("output_configuration"):
                info_msg  = "Using the old way to create the gid-output, this was removed!\n"
                info_msg += "Refer to \"https://github.com/KratosMultiphysics/Kratos/wiki/Common-"
                info_msg += "Python-Interface-of-Applications-for-Users#analysisstage-usage\" "
                info_msg += "for a description of the new format"
                raise Exception("StructuralMechanicsAnalysis: " + info_msg)
        else:
            raise NameError("wrong parameter name")

        return list_of_processes

    def _GetSimulationName(self):
        return "::[KSM Simulation]:: "

if __name__ == "__main__":
    from sys import argv

    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default ProjectParameters (read from "ProjectParameters.json"):\n'
        err_msg += '    "python3 structural_mechanics_analysis.py"\n'
        err_msg += '- With custom ProjectParameters:\n'
        err_msg += '    "python3 structural_mechanics_analysis.py CustomProjectParameters.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        project_parameters_file_name = argv[1]
    else: # using default name
        project_parameters_file_name = "ProjectParameters.json"

    with open(project_parameters_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = StructuralMechanicsAdjointDynamicAnalysis(model, parameters)
    simulation.Run()
