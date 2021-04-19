# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.kratos_utilities as KratosUtilities

# Import base class file
from KratosMultiphysics.FluidDynamicsApplication.fluid_solver import FluidSolver

class StabilizedFormulation(object):
    """Helper class to define stabilization-dependent parameters."""
    def __init__(self,settings):
        self.element_name = None
        self.element_integrates_in_time = False
        self.element_has_nodal_properties = False
        self.historical_nodal_properties_variables_list = []
        self.non_historical_nodal_properties_variables_list = []
        self.process_data = {}

        #TODO: Keep this until the MonolithicWallCondition is removed to ensure backwards compatibility in solvers with no defined condition_name
        self.condition_name = "MonolithicWallCondition"

        if settings.Has("element_type"):
            formulation = settings["element_type"].GetString()
            if formulation == "vms":
                self._SetUpClassicVMS(settings)
            elif formulation == "qsvms":
                self._SetUpQSVMS(settings)
            elif formulation == "dvms":
                self._SetUpDVMS(settings)
            elif formulation == "fic":
                self._SetUpFIC(settings)
            elif formulation == "symbolic":
                warn_msg  = 'Provided \'element_name\' is \'symbolic\'. This is has been renamed to \'weakly_compressible\'. Use this instead.'
                KratosMultiphysics.Logger.PrintWarning('\n\x1b[1;31mDEPRECATION-WARNING\x1b[0m', warn_msg)
                self._SetUpWeaklyCompressible(settings)
            elif formulation == "weakly_compressible":
                self._SetUpWeaklyCompressible(settings)
        else:
            print(settings)
            raise RuntimeError("Argument \'element_type\' not found in stabilization settings.")

    def SetProcessInfo(self,model_part):
        for variable,value in self.process_data.items():
            model_part.ProcessInfo[variable] = value

    def _SetUpClassicVMS(self,settings):
        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "vms",
            "use_orthogonal_subscales": false,
            "dynamic_tau": 0.01
        }""")

        default_non_newtonian_settings = KratosMultiphysics.Parameters(r"""{
            "power_law_k": 1e-6,
            "power_law_n": 1.0,
            "yield_stress": 0.0,
            "regularization_coefficient" : 100.0
        }""")

        # if non-newtonian, there are some extra options
        if settings.Has("non_newtonian_fluid_parameters"):
            self.non_newtonian_option = True
            default_settings.AddValue("non_newtonian_fluid_parameters", default_non_newtonian_settings)
            self.element_name = 'HerschelBulkleyVMS'
        else:
            self.non_newtonian_option = False
            self.element_name = 'VMS'
        self.condition_name = "MonolithicWallCondition"

        settings.ValidateAndAssignDefaults(default_settings)

        # set the nodal material properties flag
        self.element_has_nodal_properties = True
        self.historical_nodal_properties_variables_list = [KratosMultiphysics.VISCOSITY]

        # validate the non-newtonian parameters if necessary
        if self.non_newtonian_option:
            settings["non_newtonian_fluid_parameters"].ValidateAndAssignDefaults(default_non_newtonian_settings)
            self.process_data[KratosMultiphysics.POWER_LAW_K] = settings["non_newtonian_fluid_parameters"]["power_law_k"].GetDouble()
            self.process_data[KratosMultiphysics.POWER_LAW_N] = settings["non_newtonian_fluid_parameters"]["power_law_n"].GetDouble()
            self.process_data[KratosMultiphysics.YIELD_STRESS] = settings["non_newtonian_fluid_parameters"]["yield_stress"].GetDouble()
            self.process_data[KratosCFD.REGULARIZATION_COEFFICIENT] = settings["non_newtonian_fluid_parameters"]["regularization_coefficient"].GetDouble()

        self.process_data[KratosMultiphysics.DYNAMIC_TAU] = settings["dynamic_tau"].GetDouble()
        use_oss = settings["use_orthogonal_subscales"].GetBool()
        self.process_data[KratosMultiphysics.OSS_SWITCH] = int(use_oss)

    def _SetUpQSVMS(self,settings):
        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "qsvms",
            "use_orthogonal_subscales": false,
            "dynamic_tau": 0.0,
            "element_manages_time_integration": false
        }""")
        settings.ValidateAndAssignDefaults(default_settings)

        if settings["element_manages_time_integration"].GetBool() == False:
            self.element_name = "QSVMS"
            self.element_integrates_in_time = False
        else:
            self.element_name = "TimeIntegratedQSVMS"
            self.element_integrates_in_time = True
        self.condition_name = "NavierStokesWallCondition"

        self.process_data[KratosMultiphysics.DYNAMIC_TAU] = settings["dynamic_tau"].GetDouble()
        use_oss = settings["use_orthogonal_subscales"].GetBool()
        self.process_data[KratosMultiphysics.OSS_SWITCH] = int(use_oss)

    def _SetUpDVMS(self,settings):
        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "dvms",
            "dynamic_tau": 0.0,
            "use_orthogonal_subscales": false
        }""")
        settings.ValidateAndAssignDefaults(default_settings)

        self.element_name = "DVMS"
        self.condition_name = "NavierStokesWallCondition"

        self.process_data[KratosMultiphysics.DYNAMIC_TAU] = settings["dynamic_tau"].GetDouble()
        use_oss = settings["use_orthogonal_subscales"].GetBool()
        self.process_data[KratosMultiphysics.OSS_SWITCH] = int(use_oss)

    def _SetUpFIC(self,settings):
        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "fic",
            "beta": 0.8,
            "adjust_beta_dynamically": false,
            "dynamic_tau": 0.0
        }""")
        settings.ValidateAndAssignDefaults(default_settings)

        dynamic_beta = settings["adjust_beta_dynamically"].GetBool()
        if dynamic_beta:
            KratosMultiphysics.Logger.PrintWarning("NavierStokesSolverVMSMonolithic","FIC with dynamic beta not yet implemented, using provided beta as a constant value")
        else:
            self.element_name = "FIC"
        self.condition_name = "NavierStokesWallCondition"

        self.process_data[KratosCFD.FIC_BETA] = settings["beta"].GetDouble()
        self.process_data[KratosMultiphysics.OSS_SWITCH] = 0
        self.process_data[KratosMultiphysics.DYNAMIC_TAU] = settings["dynamic_tau"].GetDouble()

    def _SetUpWeaklyCompressible(self,settings):
        #TODO: Remove this after deprecation period is over
        if (settings["element_type"].GetString() == "symbolic"):
            settings["element_type"].SetString("weakly_compressible")

        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "weakly_compressible",
            "dynamic_tau": 1.0,
            "sound_velocity": 1.0e+5
        }""")
        settings.ValidateAndAssignDefaults(default_settings)

        self.element_name = "WeaklyCompressibleNavierStokes"
        self.condition_name = "NavierStokesWallCondition"
        self.element_integrates_in_time = True

        # set the nodal material properties flag
        self.element_has_nodal_properties = True
        self.non_historical_nodal_properties_variables_list = [KratosMultiphysics.SOUND_VELOCITY]

        self.process_data[KratosMultiphysics.DYNAMIC_TAU] = settings["dynamic_tau"].GetDouble()
        #TODO: Remove SOUND_VELOCITY from ProcessInfo. Should be obtained from the properties.
        #self.process_data[KratosMultiphysics.SOUND_VELOCITY] = settings["sound_velocity"].GetDouble()

def CreateSolver(model, custom_settings):
    return NavierStokesSolverLevelSet(model, custom_settings)

class NavierStokesSolverLevelSet(FluidSolver):

    @classmethod
    def GetDefaultParameters(cls):

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "navier_stokes_solver_levelset",
            "model_part_name": "FluidModelPart",
            "domain_size": -1,
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name",
                "reorder": false
            },
            "material_import_settings": {
                "materials_filename": ""
            },
            "formulation": {
                "element_type": "vms"
            },
            "maximum_iterations": 10,
            "echo_level": 0,
            "consider_periodic_conditions": false,
            "compute_reactions": false,
            "analysis_type": "non_linear",
            "reform_dofs_at_each_step": true,
            "relative_velocity_tolerance": 1e-3,
            "absolute_velocity_tolerance": 1e-5,
            "relative_pressure_tolerance": 1e-3,
            "absolute_pressure_tolerance": 1e-5,
            "linear_solver_settings"        : {
                "solver_type" : "amgcl"
            },
            "volume_model_part_name" : "volume_model_part",
            "skin_parts": [""],
            "assign_neighbour_elements_to_conditions": true,
            "no_skin_parts":[""],
            "time_stepping"                : {
                "automatic_time_step" : false,
                "CFL_number"          : 1,
                "minimum_delta_time"  : 1e-4,
                "maximum_delta_time"  : 0.01,
                "time_step"           : 0.0
            },
            "bfecc_convection" : false,
            "levelset_convection_settings": {
                "levelset_variable_name" : "DISTANCE",
                "levelset_convection_variable_name" : "VELOCITY",
                "levelset_gradient_variable_name" : "DISTANCE_GRADIENT",
                "max_CFL" : 1.0,
                "max_substeps" : 0,
                "levelset_splitting" : false,
                "eulerian_error_compensation" : false,
                "cross_wind_stabilization_factor" : 0.7
            },
            "distance_reinitialization" : "variational",
            "time_scheme":"bossak",
            "alpha":-0.3,
            "velocity_relaxation":0.9,
            "pressure_relaxation":0.9,
            "move_mesh_strategy": 0,
            "periodic": "periodic",
            "move_mesh_flag": false,
            "levelset_convection_settings" : {}
        }""")

        default_settings.AddMissingParameters(super(NavierStokesSolverLevelSet, cls).GetDefaultParameters())
        return default_settings

    def _BackwardsCompatibilityHelper(self,settings):
        ## Backwards compatibility -- deprecation warnings
        if settings.Has("stabilization"):
            msg  = "Input JSON data contains deprecated setting \'stabilization\'.\n"
            msg += "Please rename it to \'formulation\' (and rename \'stabilization/formulation\' to \'formulation/element_type\' if it exists).\n"
            KratosMultiphysics.Logger.PrintWarning("NavierStokesVMSMonolithicSolver",msg)
            settings.AddValue("formulation", settings["stabilization"])
            settings.RemoveValue("stabilization")
            settings["formulation"].AddValue("element_type", settings["formulation"]["formulation"])
            settings["formulation"].RemoveValue("formulation")

        if settings.Has("oss_switch"):
            msg  = "Input JSON data contains deprecated setting \'oss_switch\' (int).\n"
            msg += "Please define \'formulation/element_type\' (set it to \'vms\')\n"
            msg += "and set \'formulation/use_orthogonal_subscales\' (bool) instead."
            KratosMultiphysics.Logger.PrintWarning("NavierStokesVMSMonolithicSolver",msg)
            if not settings.Has("formulation"):
                settings.AddValue("formulation",KratosMultiphysics.Parameters(r'{"element_type":"vms"}'))
            settings["formulation"].AddEmptyValue("use_orthogonal_subscales")
            settings["formulation"]["use_orthogonal_subscales"].SetBool(bool(settings["oss_switch"].GetInt()))
            settings.RemoveValue("oss_switch")

        if settings.Has("dynamic_tau"):
            msg  = "Input JSON data contains deprecated setting \'dynamic_tau\' (float).\n"
            msg += "Please define \'formulation/element_type\' (set it to \'vms\') and \n"
            msg += "set \'formulation/dynamic_tau\' (float) instead."
            KratosMultiphysics.Logger.PrintWarning("NavierStokesVMSMonolithicSolver",msg)
            if not settings.Has("formulation"):
                settings.AddValue("formulation",KratosMultiphysics.Parameters(r'{"element_type":"vms"}'))
            settings["formulation"].AddEmptyValue("dynamic_tau")
            settings["formulation"]["dynamic_tau"].SetDouble(settings["dynamic_tau"].GetDouble())
            settings.RemoveValue("dynamic_tau")

        if settings.Has("turbulence_model") and settings["turbulence_model"].IsString():
            if settings["turbulence_model"].GetString().lower()!="none":
                msg = "Ignoring deprecated \"turbulence_model\" (string) setting."
                KratosMultiphysics.Logger.PrintWarning("NavierStokesVMSMonolithicSolver",msg)
            settings.RemoveValue("turbulence_model")

        return settings


    def __init__(self, model, custom_settings):
        self._validate_settings_in_baseclass=True # To be removed eventually
        custom_settings = self._BackwardsCompatibilityHelper(custom_settings)
        super(NavierStokesSolverLevelSet,self).__init__(model,custom_settings)

        # Set up the auxiliary class with the formulation settings
        self._SetFormulation()

        # Update the default buffer size according to the selected time scheme
        self._SetTimeSchemeBufferSize()

        self._bfecc_convection = self.settings["bfecc_convection"].GetBool()

        self._levelset_variable = KratosMultiphysics.KratosGlobals.GetVariable(self.settings["levelset_convection_settings"]["levelset_variable_name"].GetString())
        self._levelset_gradient_variable = KratosMultiphysics.KratosGlobals.GetVariable(self.settings["levelset_convection_settings"]["levelset_gradient_variable_name"].GetString())
        self._levelset_convection_variable = KratosMultiphysics.KratosGlobals.GetVariable(self.settings["levelset_convection_settings"]["levelset_convection_variable_name"].GetString())

        # this is used to perform one step BFECC for the eulerian LS convection
        self._eulerian_error_compensation = self.settings["levelset_convection_settings"]["eulerian_error_compensation"].GetBool()

        # this is used to identify the splitting of LS convection (Strang splitting idea)
        self._levelset_splitting = self.settings["levelset_convection_settings"]["levelset_splitting"].GetBool()
        self._levelset_dt_factor = 0.5 if self._levelset_splitting else 1.0
        if self._levelset_splitting:
            self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME_FACTOR, self._levelset_dt_factor)

        self._reinitialization_type = self.settings["distance_reinitialization"].GetString()
        
        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction of NavierStokesSolverLevelSet finished.")

    def AddVariables(self):
        ## Add base class variables
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_STRUCTURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BODY_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ADVPROJ)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DIVPROJ)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_WATER_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FLAG_VARIABLE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.EXTERNAL_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.Y_WALL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.Q_VALUE)

        # Adding variables required for the nodal material properties
        if self.element_has_nodal_properties:
            for variable in self.historical_nodal_properties_variables_list:
                self.main_model_part.AddNodalSolutionStepVariable(variable)

        # Adding variables required for the periodic conditions
        if self.settings["consider_periodic_conditions"].GetBool() == True:
            self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.PATCH_INDEX)

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Fluid solver variables added correctly.")

    def _GetLevelSetConvectionProcess(self):
        if not hasattr(self, '_level_set_convection_process'):
            self._level_set_convection_process = self._CreateLevelSetConvectionProcess()
        return self._level_set_convection_process

    def Initialize(self):
    	# My shock process
        settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "FluidModelPart",
            "calculate_nodal_area_at_each_step" : false,
            "shock_sensor" : true,
            "shear_sensor" : true,
            "thermal_sensor" : false,
            "thermally_coupled_formulation" : false
        }""")

        self.shock_process = KratosCFD.ShockCapturingProcess(self.model, settings)
        self.shock_process.Check()
        self.shock_process.ExecuteInitialize()
        # If the solver requires an instance of the stabilized formulation class, set the process info variables
        if hasattr(self, 'formulation'):
            self.formulation.SetProcessInfo(self.GetComputingModelPart())

        # Construct and initialize the solution strategy
        solution_strategy = self._GetSolutionStrategy()
        solution_strategy.SetEchoLevel(self.settings["echo_level"].GetInt())
        solution_strategy.Initialize()

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Solver initialization finished.")

    def InitializeSolutionStep(self):
        if self._TimeBufferIsInitialized():
            # If required, compute the BDF coefficients
            if hasattr(self, 'time_discretization'):
                (self.time_discretization).ComputeAndSaveBDFCoefficients(self.GetComputingModelPart().ProcessInfo)

            # Perform the level-set convection according to the previous step velocity
            self.__PerformLevelSetConvection()

            KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Level-set convection is performed.")

            # Recompute the distance field according to the new level-set position
            if (self._reinitialization_type == "variational"):
                self._GetDistanceReinitializationProcess().Execute()
            elif (self._reinitialization_type == "parallel"):
                adjusting_parameter = 0.05
                layers = int(adjusting_parameter*self.main_model_part.NumberOfElements()) # this parameter is essential
                max_distance = 1.0 # use this parameter to define the redistancing range
                # if using CalculateInterfacePreservingDistances(), the initial interface is preserved
                self._GetDistanceReinitializationProcess().CalculateDistances(
                    self.main_model_part,
                    self._levelset_variable,
                    KratosMultiphysics.NODAL_AREA,
                    layers,
                    max_distance,
                    self._GetDistanceReinitializationProcess().CALCULATE_EXACT_DISTANCES_TO_PLANE) #NOT_CALCULATE_EXACT_DISTANCES_TO_PLANE)

            if (self._reinitialization_type != "none"):
                KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Redistancing process is finished.")

            # distance gradient is called again to comply with the smoothed/modified DISTANCE
            self._GetDistanceGradientProcess().Execute()

            # Update the DENSITY and DYNAMIC_VISCOSITY values according to the new level-set
            self._SetNodalProperties()

            # Perform the solver InitializeSolutionStep
            self._GetSolutionStrategy().InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Mass and momentum conservation equations are solved.")
        self.shock_process.ExecuteFinalizeSolutionStep()
        if (self._levelset_splitting):
            # Perform the level-set convection to complete the solution step
            self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME_FACTOR, (1.0 - self._levelset_dt_factor))

            if self._bfecc_convection:
                self._GetLevelSetConvectionProcess().CopyScalarVarToPreviousTimeStep(
                    self.main_model_part,
                    self._levelset_variable)

            self.__PerformLevelSetConvection()


    def _SetFormulation(self):
        self.formulation = StabilizedFormulation(self.settings["formulation"])
        self.element_name = self.formulation.element_name
        self.condition_name = self.formulation.condition_name
        self.element_integrates_in_time = self.formulation.element_integrates_in_time
        self.element_has_nodal_properties = self.formulation.element_has_nodal_properties
        self.historical_nodal_properties_variables_list = self.formulation.historical_nodal_properties_variables_list
        self.non_historical_nodal_properties_variables_list = self.formulation.non_historical_nodal_properties_variables_list

    def _SetTimeSchemeBufferSize(self):
        scheme_type = self.settings["time_scheme"].GetString()
        if scheme_type == "bossak":
            self.min_buffer_size = 2
        elif scheme_type == "bdf2":
            self.min_buffer_size = 3
        elif scheme_type == "steady":
            self.min_buffer_size = 1
            self._SetUpSteadySimulation()
        else:
            msg  = "Unknown time_scheme option found in project parameters:\n"
            msg += "\"" + scheme_type + "\"\n"
            msg += "Accepted values are \"bossak\", \"bdf2\" or \"steady\".\n"
            raise Exception(msg)

    def _SetUpSteadySimulation(self):
        '''Overwrite time stepping parameters so that they do not interfere with steady state simulations.'''
        self.settings["time_stepping"]["automatic_time_step"].SetBool(False)
        if self.settings["formulation"].Has("dynamic_tau"):
            self.settings["formulation"]["dynamic_tau"].SetDouble(0.0)

    def __PerformLevelSetConvection(self):
        if self._bfecc_convection:
            self._GetLevelsetGradientProcess().Execute() #Level-set gradient is needed for the limiter
            self._GetLevelSetConvectionProcess().BFECCconvect(
                self.main_model_part,
                self._levelset_variable,
                self._levelset_convection_variable,
                self.settings["bfecc_number_substeps"].GetInt())
        else:
            if (self._eulerian_error_compensation):
                self._GetLevelsetGradientProcess().Execute() #Level-set gradient is needed for the limiter
            self._GetLevelSetConvectionProcess().Execute()


    def _SetNodalProperties(self):
        set_viscosity = KratosMultiphysics.VISCOSITY in self.historical_nodal_properties_variables_list
        set_sound_velocity = KratosMultiphysics.SOUND_VELOCITY in self.non_historical_nodal_properties_variables_list
        # Get density and dynamic viscostity from the properties of the first element
        for el in self.main_model_part.Elements:
            # Get DENSITY from properties
            # if set_density:
            #     rho = el.Properties.GetValue(KratosMultiphysics.DENSITY)
            #     if rho <= 0.0:
            #         raise Exception("DENSITY set to {0} in Properties {1}, positive number expected.".format(rho,el.Properties.Id))
            # Get DYNAMIC_VISCOSITY from properties and calculate the kinematic one (VISCOSITY)
            if set_viscosity:
                dyn_viscosity = el.Properties.GetValue(KratosMultiphysics.DYNAMIC_VISCOSITY)
                if dyn_viscosity <= 0.0:
                    raise Exception("DYNAMIC_VISCOSITY set to {0} in Properties {1}, positive number expected.".format(dyn_viscosity,el.Properties.Id))
                kin_viscosity = dyn_viscosity / rho
            # Get SOUND_VELOCITY
            # if set_sound_velocity:
            #     if el.Properties.Has(KratosMultiphysics.SOUND_VELOCITY):
            #         sound_velocity = el.Properties.GetValue(KratosMultiphysics.SOUND_VELOCITY)
            #     else:
            #         sound_velocity = 1.0e+12 # Default sound velocity value
            #         KratosMultiphysics.Logger.PrintWarning('No \'SOUND_VELOCITY\' value fod in Properties {0}. Setting default value {1}'.format(el.Properties.Id, sound_velocity))
            #     if sound_velocity <= 0.0:
            #         raise Exception("SOUND_VELOCITY set to {0} in Properties {1}, positive number expected.".format(sound_velocity, el.Properties.Id))
            break
        else:
            raise Exception("No fluid elements found in the main model part.")

        # Transfer the obtained properties to the nodes
        # if set_density:
        #     KratosMultiphysics.VariableUtils().SetVariable(KratosMultiphysics.DENSITY, rho, self.main_model_part.Nodes)
        if set_viscosity:
            KratosMultiphysics.VariableUtils().SetVariable(KratosMultiphysics.VISCOSITY, kin_viscosity, self.main_model_part.Nodes)
<<<<<<< HEAD
        # if set_sound_velocity:
        #     KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.SOUND_VELOCITY, sound_velocity, self.main_model_part.Nodes)
    
=======
        #if set_sound_velocity:
            #KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.SOUND_VELOCITY, sound_velocity, self.main_model_part.Nodes)

>>>>>>> origin/level-set-area-branch
    def _CreateLevelSetConvectionProcess(self):
        # Construct the level set convection process
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        computing_model_part = self.GetComputingModelPart()
        if self._bfecc_convection:
            if have_conv_diff:
                if domain_size == 2:
                    locator = KratosMultiphysics.BinBasedFastPointLocator2D(computing_model_part)
                    locator.UpdateSearchDatabase()
                    level_set_convection_process = KratosConvDiff.BFECCConvection2D(locator, self._levelset_splitting, True)
                else:
                    locator = KratosMultiphysics.BinBasedFastPointLocator3D(computing_model_part)
                    locator.UpdateSearchDatabase()
                    level_set_convection_process = KratosConvDiff.BFECCConvection3D(locator, self._levelset_splitting, True)
            else:
                raise Exception("The BFECC level set convection requires the Kratos ConvectionDiffusionApplication compilation.")
        else:
            linear_solver = self._GetLinearSolver()

            levelset_convection_settings = self.settings["levelset_convection_settings"]
            if domain_size == 2:
                level_set_convection_process = KratosMultiphysics.LevelSetConvectionProcess2D(
                    computing_model_part,
                    linear_solver,
                    levelset_convection_settings)
            else:
                level_set_convection_process = KratosMultiphysics.LevelSetConvectionProcess3D(
                    computing_model_part,
                    linear_solver,
                    levelset_convection_settings)

        return level_set_convection_process

    def _CreateLevelsetGradientProcess(self):
        levelset_gradient_process = KratosMultiphysics.ComputeNonHistoricalNodalGradientProcess(
                self.main_model_part,
                self._levelset_variable,
                self._levelset_gradient_variable,
                KratosMultiphysics.NODAL_AREA)

        return levelset_gradient_process

    def _GetLevelsetGradientProcess(self):
        if not hasattr(self, '_levelset_gradient_process'):
            self._levelset_gradient_process = self._CreateLevelsetGradientProcess()
        return self._levelset_gradient_process

    def _GetDistanceReinitializationProcess(self):
        if not hasattr(self, '_distance_reinitialization_process'):
            self._distance_reinitialization_process = self._CreateDistanceReinitializationProcess()
        return self._distance_reinitialization_process

    def _CreateDistanceReinitializationProcess(self):
        # Construct the variational distance calculation process
        if (self._reinitialization_type == "variational"):
            maximum_iterations = 2 #TODO: Make this user-definable
            linear_solver = self._GetLinearSolver()
            computing_model_part = self.GetComputingModelPart()
            if self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
                distance_reinitialization_process = KratosMultiphysics.VariationalDistanceCalculationProcess2D(
                    computing_model_part,
                    linear_solver,
                    maximum_iterations,
                    KratosMultiphysics.VariationalDistanceCalculationProcess2D.CALCULATE_EXACT_DISTANCES_TO_PLANE)
            else:
                distance_reinitialization_process = KratosMultiphysics.VariationalDistanceCalculationProcess3D(
                    computing_model_part,
                    linear_solver,
                    maximum_iterations,
                    KratosMultiphysics.VariationalDistanceCalculationProcess3D.CALCULATE_EXACT_DISTANCES_TO_PLANE)

        elif (self._reinitialization_type == "parallel"):
            if self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
                locator = KratosMultiphysics.BinBasedFastPointLocator2D(self.main_model_part)
                locator.UpdateSearchDatabase()
                distance_reinitialization_process = KratosMultiphysics.ParallelDistanceCalculator2D()
            else:
                locator = KratosMultiphysics.BinBasedFastPointLocator3D(self.main_model_part)
                locator.UpdateSearchDatabase()
                distance_reinitialization_process = KratosMultiphysics.ParallelDistanceCalculator3D()
        elif (self._reinitialization_type == "none"):
                KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Redistancing is turned off.")
        else:
            raise Exception("Please use a valid distance reinitialization type or set it as \'none\'. Valid types are: \'variational\' and \'parallel\'.")

        return distance_reinitialization_process

    def _GetDistanceGradientProcess(self):
        if not hasattr(self, '_distance_gradient_process'):
            self._distance_gradient_process = self._CreateDistanceGradientProcess()
        return self._distance_gradient_process

    def _CreateDistanceGradientProcess(self):
        distance_gradient_process = KratosMultiphysics.ComputeNodalGradientProcess(
                self.main_model_part,
                self._levelset_variable,
                self._levelset_gradient_variable,
                KratosMultiphysics.NODAL_AREA)

        return distance_gradient_process