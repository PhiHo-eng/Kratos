# Importing Kratos
import KratosMultiphysics
from KratosMultiphysics.process_factory import KratosProcessFactory
from KratosMultiphysics.kratos_utilities import IssueDeprecationWarning

class AnalysisStage(object):
    """The base class for the AnalysisStage-classes in the applications
    Changes to this BaseClass have to be discussed first!
    """
    def __init__(self, model, project_parameters):
        """The constructor of the AnalysisStage-Object.

        It is intended to be called from the constructor
        of deriving classes:
        super(DerivedAnalysis, self).__init__(project_parameters)

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- The Model to be used
        project_parameters -- The ProjectParameters used
        """
        if not isinstance(model, KratosMultiphysics.Model):
            raise Exception("Input is expected to be provided as a Kratos Model object")

        if not isinstance(project_parameters, KratosMultiphysics.Parameters):
            raise Exception("Input is expected to be provided as a Kratos Parameters object")

        self.model = model
        self.project_parameters = project_parameters

        # We validate reinitilize settings
        reinitialize_settings = self._GetReInitializeRequired()
        if not self.project_parameters.Has("reinitialize_settings"):
            self.project_parameters.AddEmptyValue("reinitialize_settings")
        self.project_parameters["reinitialize_settings"].ValidateAndAssignDefaults(reinitialize_settings)

        ## Get echo level and parallel type
        self.echo_level = self.project_parameters["problem_data"]["echo_level"].GetInt()
        self.parallel_type = self.project_parameters["problem_data"]["parallel_type"].GetString()
        is_distributed_run = KratosMultiphysics.IsDistributedRun()

        if self.parallel_type == "OpenMP" and is_distributed_run:
            KratosMultiphysics.Logger.PrintWarning("Parallel Type", '"OpenMP" is specified as "parallel_type", but Kratos is running distributed!')
        if self.parallel_type == "MPI" and not is_distributed_run:
            KratosMultiphysics.Logger.PrintWarning("Parallel Type", '"MPI" is specified as "parallel_type", but Kratos is not running distributed!')

        self._GetSolver().AddVariables() # this creates the solver and adds the variables

    def Run(self):
        """This function executes the entire AnalysisStage
        It can be overridden by derived classes
        """
        self.Initialize()
        self.RunSolutionLoop()
        self.Finalize()

    def KeepAdvancingSolutionLoop(self):
        """This function specifies the stopping criteria for breaking the solution loop
        It can be overridden by derived classes
        """
        return self.time < self.end_time

    def RunSolutionLoop(self):
        """This function executes the solution loop of the AnalysisStage
        It can be overridden by derived classes
        """
        while self.KeepAdvancingSolutionLoop():
            self.time = self._GetSolver().AdvanceInTime(self.time)
            self.InitializeSolutionStep()
            self._GetSolver().Predict()
            is_converged = self._GetSolver().SolveSolutionStep()
            self.__CheckIfSolveSolutionStepReturnsAValue(is_converged)
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()

    def Initialize(self):
        """This function initializes the AnalysisStage
        Usage: It is designed to be called ONCE, BEFORE the execution of the solution-loop
        This function has to be implemented in deriving classes!
        """
        # Modelers:
        self._CreateModelers()
        self._ModelersSetupGeometryModel()
        self._ModelersPrepareGeometryModel()
        self._ModelersSetupModelPart()

        self._GetSolver().ImportModelPart()
        self._GetSolver().PrepareModelPart()
        self._GetSolver().AddDofs()

        self.ModifyInitialProperties()
        self.ModifyInitialGeometry()

        ##here we initialize user-provided processes
        self.__CreateListOfProcesses() # has to be done after importing and preparing the ModelPart
        for process in self._GetListOfProcesses():
            process.ExecuteInitialize()

        self._GetSolver().Initialize()
        self.Check()

        self.ModifyAfterSolverInitialize()

        for process in self._GetListOfProcesses():
            process.ExecuteBeforeSolutionLoop()

        ## Stepping and time settings
        self.end_time = self.project_parameters["problem_data"]["end_time"].GetDouble()

        if self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            self.time = self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME]
        else:
            self.time = self.project_parameters["problem_data"]["start_time"].GetDouble()

        ## If the echo level is high enough, print the complete list of settings used to run the simualtion
        if self.echo_level > 1:
            with open("ProjectParametersOutput.json", 'w') as parameter_output_file:
                parameter_output_file.write(self.project_parameters.PrettyPrintJsonString())

        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Analysis -START- ")

    def Finalize(self):
        """This function finalizes the AnalysisStage
        Usage: It is designed to be called ONCE, AFTER the execution of the solution-loop
        """
        for process in self._GetListOfProcesses():
            process.ExecuteFinalize()

        self._GetSolver().Finalize()

        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Analysis -END- ")

    def InitializeSolutionStep(self):
        """This function performs all the required operations that should be executed
        (for each step) BEFORE solving the solution step.
        """
        # We reinitialize if remeshed previously
        if self._CheckIfModelIsModified():
            self.ClearDatabase()
            self.ReInitializeSolver()

        self.PrintAnalysisStageProgressInformation()

        self.ApplyBoundaryConditions() #here the processes are called
        self.ChangeMaterialProperties() #this is normally empty
        self._GetSolver().InitializeSolutionStep()

        # We reinitialize if remeshed on the InitializeSolutionStep
        if self._CheckIfModelIsModified():
            self.ClearDatabase()
            self.ReInitializeSolver()
            self.InitializeSolutionStep()
        
    def PrintAnalysisStageProgressInformation(self):
        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "STEP: ", self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.STEP])
        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "TIME: ", self.time)

    def FinalizeSolutionStep(self):
        """This function performs all the required operations that should be executed
        (for each step) AFTER solving the solution step.
        """
        self._GetSolver().FinalizeSolutionStep()

        for process in self._GetListOfProcesses():
            process.ExecuteFinalizeSolutionStep()

    def OutputSolutionStep(self):
        """This function printed / writes output files after the solution of a step
        """
        execute_was_called = False
        for output_process in self._GetListOfOutputProcesses():
            if output_process.IsOutputStep():
                if not execute_was_called:
                    for process in self._GetListOfProcesses():
                        process.ExecuteBeforeOutputStep()
                    execute_was_called = True

                output_process.PrintOutput()

        if execute_was_called:
            for process in self._GetListOfProcesses():
                process.ExecuteAfterOutputStep()

    def Check(self):
        """This function checks the AnalysisStage
        """
        self._GetSolver().Check()
        for process in self._GetListOfProcesses():
            process.Check()

    def Clear(self):
        """This function clears the AnalysisStage
        """
        self._GetSolver().Clear()

    def ModifyInitialProperties(self):
        """this is the place to eventually modify material properties in the stage """
        pass

    def ModifyInitialGeometry(self):
        """this is the place to eventually modify geometry (for example moving nodes) in the stage """
        pass

    def ModifyAfterSolverInitialize(self):
        """this is the place to eventually do any modification that requires the solver to be initialized """
        pass

    def ApplyBoundaryConditions(self):
        """here the boundary conditions is applied, by calling the InitializeSolutionStep function of the processes"""

        for process in self._GetListOfProcesses():
            process.ExecuteInitializeSolutionStep()

        #other operations as needed

    def ChangeMaterialProperties(self):
        """this function is where the user could change material parameters as a part of the solution step """
        pass

    def ClearDatabase(self):
        """ This method clears the database in case it is necessary. For example it is used in the remeshing process when some specific modifications must be done in the database (i.e To reset a flag that affects the whole model part)
            Keyword arguments:
            self It signifies an instance of a class.
        """
        pass
    
    def ReInitializeSolver(self):
        """ This reinitializes the solver and the processes (used for example on remesh or on adaptive NR)

            Keyword arguments:
            self It signifies an instance of a class.
        """
        solver = self._GetSolver()
        processes = self._GetListOfProcesses()

        reinitialize_settings = self.project_parameters["reinitialize_settings"]

        # Clear solver
        if reinitialize_settings["clear_solver"].GetBool():
            solver.Clear()

        # Some initial calls
        if reinitialize_settings["modify_initial_properties"].GetBool():
            self.ModifyInitialProperties()
        if reinitialize_settings["modify_initial_geometries"].GetBool():
            self.ModifyInitialGeometry()

        # Processes initialization
        if reinitialize_settings["initialize_processes"].GetBool():
            for process in processes:
                process.ExecuteInitialize()

            if not reinitialize_settings["before_solution_loop_processes"].GetBool():
                KratosMultiphysics.Logger.PrintWarning(self._GetSimulationName(), "If ExecuteInitialize in processes is executed, ExecuteBeforeSolutionLoop should be called")

        # WE INITIALIZE THE SOLVER
        if reinitialize_settings["solver_initialize"].GetBool():
            solver.Initialize()
        elif reinitialize_settings["entities_initialize"].GetBool(): # Solver already initializes entities
            # Initilize elements and conditions
            model_part_names = self.model.GetModelPartNames()
            for model_part_name in model_part_names:
                KratosMultiphysics.EntitiesUtilities.InitializeAllEntities(self.model.GetModelPart(model_part_name))

        # Call check
        if reinitialize_settings["check_analysis"].GetBool():
            self.Check()

        # Modify after initialize
        if reinitialize_settings["modify_after_solver_initialize"].GetBool():
            self.ModifyAfterSolverInitialize()

        ## Processes before the loop
        if reinitialize_settings["before_solution_loop_processes"].GetBool():
            for process in processes:
                process.ExecuteBeforeSolutionLoop()

            if not reinitialize_settings["initialize_solution_step_processes"].GetBool():
                KratosMultiphysics.Logger.PrintWarning(self._GetSimulationName(), "If ExecuteBeforeSolutionLoop in processes is executed, ExecuteInitializeSolutionStep should be called")

        ## Processes of initialize the solution step
        if reinitialize_settings["initialize_solution_step_processes"].GetBool():
            for process in processes:
                process.ExecuteInitializeSolutionStep()

        # We reset the flags
        self._ResetModelIsModified()
      
    def _GetReInitializeRequired(self):
        """ This returns the initilization requirement. By default only elements and conditions are initialized

            Keyword arguments:
            self It signifies an instance of a class.
        """

        reinitialize_settings = KratosMultiphysics.Parameters("""{
            "clear_solver"                       : true,
            "modify_initial_properties"          : false,
            "modify_initial_geometries"          : false,
            "initialize_processes"               : true,
            "solver_initialize"                  : false,
            "entities_initialize"                : true,
            "check_analysis"                     : false,
            "modify_after_solver_initialize"     : false,
            "before_solution_loop_processes"     : true,
            "initialize_solution_step_processes" : true
        }""")

        return reinitialize_settings
      
    def _GetSolver(self):
        if not hasattr(self, '_solver'):
            self._solver = self._CreateSolver()
        return self._solver

    def _CreateSolver(self):
        """Create the solver
        """
        raise Exception("Creation of the solver must be implemented in the derived class.")

    ### Modelers
    def _ModelersSetupGeometryModel(self):
        # Import or generate geometry models from external input.
        for modeler in self._GetListOfModelers():
            if self.echo_level > 1:
                KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Modeler: ", str(modeler), " Setup Geometry Model started.")
            modeler.SetupGeometryModel()
            if self.echo_level > 1:
                KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Modeler: ", str(modeler), " Setup Geometry Model finished.")

    def _ModelersPrepareGeometryModel(self):
        # Prepare or update the geometry model_part.
        for modeler in self._GetListOfModelers():
            if self.echo_level > 1:
                KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Modeler: ", str(modeler), " Prepare Geometry Model started.")
            modeler.PrepareGeometryModel()
            if self.echo_level > 1:
                KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Modeler: ", str(modeler), " Prepare Geometry Model finished.")

    def _ModelersSetupModelPart(self):
        # Convert the geometry model or import analysis suitable models.
        for modeler in self._GetListOfModelers():
            if self.echo_level > 1:
                KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Modeler: ", str(modeler), " Setup ModelPart started.")
            modeler.SetupModelPart()
            if self.echo_level > 1:
                KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Modeler: ", str(modeler), " Setup ModelPart finished.")

    ### Modelers
    def _GetListOfModelers(self):
        """ This function returns the list of modelers
        """
        if not hasattr(self, '_list_of_modelers'):
            raise Exception("The list of modelers was not yet created!")
        return self._list_of_modelers

    def _CreateModelers(self):
        """ List of modelers in following format:
        "modelers" : [{
            "modeler_name" : "geometry_import":
            "parameters" : {
                "echo_level" : 0:
                // settings for this modeler
            }
        },{ ... }]
        """
        self._list_of_modelers = []

        if self.project_parameters.Has("modelers"):
            from KratosMultiphysics.modeler_factory import KratosModelerFactory
            factory = KratosModelerFactory()

            modelers_list = self.project_parameters["modelers"]
            self._list_of_modelers = factory.ConstructListOfModelers(self.model, modelers_list)

    ### Processes
    def _GetListOfProcesses(self):
        """This function returns the list of processes involved in this Analysis
        """
        if not hasattr(self, '_list_of_processes'):
            raise Exception("The list of processes was not yet created!")
        return self._list_of_processes

    def _GetListOfOutputProcesses(self):
        """This function returns the list of output processes involved in this Analysis
        """
        if not hasattr(self, '_list_of_output_processes'):
            raise Exception("The list of output-processes was not yet created!")
        return self._list_of_output_processes

    def _CreateProcesses(self, parameter_name, initialization_order):
        """Create a list of processes
        Format:
        "processes" : {
            initial_processes : [
                { proces_specific_params },
                { proces_specific_params }
            ],
            boundary_processes : [
                { proces_specific_params },
                { proces_specific_params }
            ]
        }
        The order of intialization can be specified by setting it in "initialization_order"
        if e.g. the "boundary_processes" should be constructed before the "initial_processes", then
        initialization_order should be a list containing ["boundary_processes", "initial_processes"]
        see the functions _GetOrderOfProcessesInitialization and _GetOrderOfOutputProcessesInitialization
        """
        list_of_processes = []

        factory = KratosProcessFactory(self.model)

        if self.project_parameters.Has(parameter_name):
            processes_params = self.project_parameters[parameter_name]

            # first initialize the processes that depend on the order
            for processes_names in initialization_order:
                if processes_params.Has(processes_names):
                    list_of_processes += factory.ConstructListOfProcesses(processes_params[processes_names])

            # then initialize the processes that don't depend on the order
            for name, value in processes_params.items():
                if not name in initialization_order:
                    list_of_processes += factory.ConstructListOfProcesses(value) # Does this work? or should it be processes[name]

        return list_of_processes

    def _GetOrderOfProcessesInitialization(self):
        """This function can be overridden in derived classes if the order of
        initialization for the processes matters
        """
        return []

    def _GetOrderOfOutputProcessesInitialization(self):
        """This function can be overridden in derived classes if the order of
        initialization for the output-processes matters
        """
        return []

    def _GetSimulationName(self):
        """Returns the name of the Simulation
        """
        return "Analysis"

    def __CreateListOfProcesses(self):
        """This function creates the processes and the output-processes
        """
        order_processes_initialization = self._GetOrderOfProcessesInitialization()
        self._list_of_processes        = self._CreateProcesses("processes", order_processes_initialization)
        order_processes_initialization = self._GetOrderOfOutputProcessesInitialization()
        self._list_of_output_processes = self._CreateProcesses("output_processes", order_processes_initialization)
        self._list_of_processes.extend(self._list_of_output_processes) # Adding the output processes to the regular processes

    def __CheckIfSolveSolutionStepReturnsAValue(self, is_converged):
        """In case the solver does not return the state of convergence
        (same as the SolvingStrategy does) then issue ONCE a deprecation-warning

        """
        if is_converged is None:
            if not hasattr(self, '_map_ret_val_depr_warnings'):
                self._map_ret_val_depr_warnings = []
            solver_class_name = self._GetSolver().__class__.__name__
            # used to only print the deprecation-warning once
            if not solver_class_name in self._map_ret_val_depr_warnings:
                self._map_ret_val_depr_warnings.append(solver_class_name)
                warn_msg  = 'Solver "{}" does not return '.format(solver_class_name)
                warn_msg += 'the state of convergence from "SolveSolutionStep"'
                IssueDeprecationWarning("AnalysisStage", warn_msg)
                
    def _CheckIfModelIsModified(self):
        """ This checks the flag MODIFIED in all the modelparts belonging to the analysis model. Returns true if at least one of the model parts is modified
            Keyword arguments:
            self It signifies an instance of a class.
        """
        modified_model = False
        if hasattr(self, 'model'):
            model_part_names = self.model.GetModelPartNames()
            for model_part_name in model_part_names:
                if self.model.GetModelPart(model_part_name).Is(KratosMultiphysics.MODIFIED):
                    modified_model = True
                    break

        return modified_model

    def _SetModelIsModified(self):
        """ This sets the flag MODIFIED in all the modelparts belonging to the analysis model
            Keyword arguments:
            self It signifies an instance of a class.
        """
        if hasattr(self, 'model'):
            model_part_names = self.model.GetModelPartNames()
            for model_part_name in model_part_names:
                self.model.GetModelPart(model_part_name).Set(KratosMultiphysics.MODIFIED, True)

    def _ResetModelIsModified(self):
        """ This resets the flag MODIFIED in all the modelparts belonging to the analysis model
            Keyword arguments:
            self It signifies an instance of a class.
        """
        if hasattr(self, 'model'):
            model_part_names = self.model.GetModelPartNames()
            for model_part_name in model_part_names:
                self.model.GetModelPart(model_part_name).Reset(KratosMultiphysics.MODIFIED)
