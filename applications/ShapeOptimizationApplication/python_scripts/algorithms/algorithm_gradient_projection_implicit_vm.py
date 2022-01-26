# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Reza Najian Asl
#                   
#
# ==============================================================================

# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Kratos Core and Apps
import KratosMultiphysics as KM
import KratosMultiphysics.ShapeOptimizationApplication as KSO
from KratosMultiphysics.LinearSolversApplication import dense_linear_solver_factory

# Additional imports
from KratosMultiphysics.ShapeOptimizationApplication.algorithms.algorithm_base import OptimizationAlgorithm
from KratosMultiphysics.ShapeOptimizationApplication import mapper_factory
from KratosMultiphysics.ShapeOptimizationApplication.loggers import data_logger_factory
from KratosMultiphysics.ShapeOptimizationApplication.utilities.custom_timer import Timer
from KratosMultiphysics.ShapeOptimizationApplication.utilities.custom_variable_utilities import WriteDictionaryDataOnNodalVariable
import math

# ==============================================================================
class AlgorithmGradientProjectionImplicitVM(OptimizationAlgorithm):
    # --------------------------------------------------------------------------
    def __init__(self, optimization_settings, analyzer, communicator, model_part_controller):
        default_algorithm_settings = KM.Parameters("""
        {
            "name"                    : "gradient_projection_implicit_vm",
            "max_correction_share"    : 0.75,
            "max_iterations"          : 100,
            "relative_tolerance"      : 1e-3,
            "line_search" : {
                "line_search_type"           : "manual_stepping",
                "normalize_search_direction" : true,
                "step_size"                  : 1.0
            }
        }""")
        self.algorithm_settings =  optimization_settings["optimization_algorithm"]
        self.algorithm_settings.RecursivelyValidateAndAssignDefaults(default_algorithm_settings)

        self.optimization_settings = optimization_settings
        self.mapper_settings = optimization_settings["design_variables"]["filter"]

        self.analyzer = analyzer
        self.communicator = communicator
        self.model_part_controller = model_part_controller

        self.design_surface = None
        self.mapper = None
        self.data_logger = None
        self.optimization_utilities = None

        self.objectives = optimization_settings["objectives"]
        self.constraints = optimization_settings["constraints"]
        self.constraint_gradient_variables = {}
        for itr, constraint in enumerate(self.constraints):
            self.constraint_gradient_variables.update({
                constraint["identifier"].GetString() : {
                    "gradient": KM.KratosGlobals.GetVariable("DC"+str(itr+1)+"DX"),
                    "mapped_gradient": KM.KratosGlobals.GetVariable("DC"+str(itr+1)+"DX_MAPPED")
                }
            })
        self.max_correction_share = self.algorithm_settings["max_correction_share"].GetDouble()

        self.step_size = self.algorithm_settings["line_search"]["step_size"].GetDouble()
        self.max_iterations = self.algorithm_settings["max_iterations"].GetInt() + 1
        self.relative_tolerance = self.algorithm_settings["relative_tolerance"].GetDouble()

        self.optimization_model_part = model_part_controller.GetOptimizationModelPart()
        self.optimization_model_part.AddNodalSolutionStepVariable(KSO.SEARCH_DIRECTION)
        self.optimization_model_part.AddNodalSolutionStepVariable(KSO.CORRECTION)
        self.optimization_model_part.AddNodalSolutionStepVariable(KSO.HELMHOLTZ_VARS)
        self.optimization_model_part.AddNodalSolutionStepVariable(KSO.HELMHOLTZ_SOURCE)
        self.optimization_model_part.AddNodalSolutionStepVariable(KSO.CONTROL_POINT) 
        self.optimization_model_part.AddNodalSolutionStepVariable(KSO.SHAPE)       

    # --------------------------------------------------------------------------
    def CheckApplicability(self):
        if self.objectives.size() > 1:
            raise RuntimeError("Gradient projection algorithm only supports one objective function!")
        if self.constraints.size() == 0:
            raise RuntimeError("Gradient projection algorithm requires definition of at least one constraint!")

    # --------------------------------------------------------------------------
    def InitializeOptimizationLoop(self):
        self.model_part_controller.Initialize()

        self.design_surface = self.model_part_controller.GetDesignSurface()

        self.mapper = mapper_factory.CreateMapper(self.optimization_model_part, self.design_surface, self.mapper_settings)

        # here we add DOFs needed by implicit vertex-morphing.
        if(self.mapper.implicit_settings["only_design_surface_parameterization"].GetBool()):
            for node in self.design_surface.Nodes:
                node.AddDof(KSO.HELMHOLTZ_VARS_X)
                node.AddDof(KSO.HELMHOLTZ_VARS_Y)
                node.AddDof(KSO.HELMHOLTZ_VARS_Z)
            for key, value in self.model_part_controller.damping_regions.items():
                for node in value.Nodes:
                    node.AddDof(KSO.HELMHOLTZ_VARS_X)
                    node.AddDof(KSO.HELMHOLTZ_VARS_Y)
                    node.AddDof(KSO.HELMHOLTZ_VARS_Z)                
        else:
            for node in self.optimization_model_part.Nodes:
                node.AddDof(KSO.HELMHOLTZ_VARS_X)
                node.AddDof(KSO.HELMHOLTZ_VARS_Y)
                node.AddDof(KSO.HELMHOLTZ_VARS_Z)            

        # here we fixed the damping/fixed regions        
        index=0
        for key, value in self.model_part_controller.damping_regions.items():
            for node in value.Nodes:
                if(self.model_part_controller.model_settings["damping"]["damping_regions"][index]["damp_X"].GetBool()):
                    node.Fix(KSO.HELMHOLTZ_VARS_X)
                if(self.model_part_controller.model_settings["damping"]["damping_regions"][index]["damp_Y"].GetBool()):
                    node.Fix(KSO.HELMHOLTZ_VARS_Y)            
                if(self.model_part_controller.model_settings["damping"]["damping_regions"][index]["damp_Z"].GetBool()):
                    node.Fix(KSO.HELMHOLTZ_VARS_Z)   
            index += 1
 
        # here we save the initial coords
        self.init_coords = []
        for node in self.optimization_model_part.Nodes:
            self.init_coords.append(node.X0)
            self.init_coords.append(node.Y0)
            self.init_coords.append(node.Z0)             

        self.analyzer.InitializeBeforeOptimizationLoop()

        self.model_part_controller.ComputeUnitSurfaceNormals()
        
        self.mapper.Initialize()

        self.data_logger = data_logger_factory.CreateDataLogger(self.model_part_controller, self.communicator, self.optimization_settings)
        self.data_logger.InitializeDataLogging()

        self.optimization_utilities = KSO.OptimizationUtilities

    # --------------------------------------------------------------------------
    def RunOptimizationLoop(self):
        timer = Timer()
        timer.StartTimer()

        for self.optimization_iteration in range(1,self.max_iterations):
            KM.Logger.Print("")
            KM.Logger.Print("===============================================================================")
            KM.Logger.PrintInfo("ShapeOpt", timer.GetTimeStamp(), ": Starting optimization iteration ", self.optimization_iteration)
            KM.Logger.Print("===============================================================================\n")

            timer.StartNewLap()

            self.__initializeNewShape()

            self.__analyzeShape()

            self.__computeShapeUpdate()

            self.__logCurrentOptimizationStep()

            KM.Logger.Print("")
            KM.Logger.PrintInfo("ShapeOpt", "Time needed for current optimization step = ", timer.GetLapTime(), "s")
            KM.Logger.PrintInfo("ShapeOpt", "Time needed for total optimization so far = ", timer.GetTotalTime(), "s")

            if self.__isAlgorithmConverged():
                break
            else:
                self.__determineAbsoluteChanges()

    # --------------------------------------------------------------------------
    def FinalizeOptimizationLoop(self):
        self.data_logger.FinalizeDataLogging()
        self.analyzer.FinalizeAfterOptimizationLoop()

    # --------------------------------------------------------------------------
    def __initializeNewShape(self):
        self.model_part_controller.UpdateTimeStep(self.optimization_iteration)

        # update the control points
        for node in self.optimization_model_part.Nodes:
            node_control_point = node.GetSolutionStepValue(KSO.CONTROL_POINT)
            node_control_point_update = node.GetSolutionStepValue(KSO.CONTROL_POINT_UPDATE)
            node_control_point += node_control_point_update
            node.SetSolutionStepValue(KSO.CONTROL_POINT,node_control_point)        
        # then update the shape  

        if(self.mapper.implicit_settings["only_design_surface_parameterization"].GetBool()):              
            self.model_part_controller.UpdateMeshAccordingInputVariable(KSO.SHAPE_UPDATE)
            self.model_part_controller.SetReferenceMeshToMesh()
        else:
            for node in self.optimization_model_part.Nodes:
                node_shape_update = node.GetSolutionStepValue(KSO.SHAPE_UPDATE)  
                node.X += node_shape_update[0]
                node.X0 += node_shape_update[0]
                node.Y += node_shape_update[1]
                node.Y0 += node_shape_update[1]
                node.Z += node_shape_update[2]
                node.Z0 += node_shape_update[2]  
        
        self.model_part_controller.ComputeUnitSurfaceNormals()       

    # --------------------------------------------------------------------------
    def __analyzeShape(self):
        self.communicator.initializeCommunication()
        self.communicator.requestValueOf(self.objectives[0]["identifier"].GetString())
        self.communicator.requestGradientOf(self.objectives[0]["identifier"].GetString())

        for constraint in self.constraints:
            con_id =  constraint["identifier"].GetString()
            self.communicator.requestValueOf(con_id)
            self.communicator.requestGradientOf(con_id)

        self.analyzer.AnalyzeDesignAndReportToCommunicator(self.optimization_model_part, self.optimization_iteration, self.communicator)

        # compute normals only if required
        surface_normals_required = self.objectives[0]["project_gradient_on_surface_normals"].GetBool()
        for constraint in self.constraints:
            if constraint["project_gradient_on_surface_normals"].GetBool():
                surface_normals_required = True

        if surface_normals_required:
            self.model_part_controller.ComputeUnitSurfaceNormals()

        # project and damp objective gradients
        objGradientDict = self.communicator.getStandardizedGradient(self.objectives[0]["identifier"].GetString())
        WriteDictionaryDataOnNodalVariable(objGradientDict, self.optimization_model_part, KSO.DF1DX)

        if self.objectives[0]["project_gradient_on_surface_normals"].GetBool():
            self.model_part_controller.ProjectNodalVariableOnUnitSurfaceNormals(KSO.DF1DX)

        if self.objectives[0]["deintegrate_area"].GetBool():
            self.model_part_controller.ComputeUnitSurfaceNormals()
            self.model_part_controller.AreaDeintegrateNodalVariable(KSO.DF1DX) 

        # project and damp constraint gradients
        for constraint in self.constraints:
            con_id = constraint["identifier"].GetString()
            conGradientDict = self.communicator.getStandardizedGradient(con_id)
            gradient_variable = self.constraint_gradient_variables[con_id]["gradient"]
            WriteDictionaryDataOnNodalVariable(conGradientDict, self.optimization_model_part, gradient_variable)

            if constraint["project_gradient_on_surface_normals"].GetBool():
                self.model_part_controller.ProjectNodalVariableOnUnitSurfaceNormals(gradient_variable)
            if constraint["deintegrate_area"].GetBool():
                self.model_part_controller.AreaDeintegrateNodalVariable(gradient_variable)                

    # --------------------------------------------------------------------------
    def __computeShapeUpdate(self):

        self.mapper.Update()

        if(self.mapper.implicit_settings["formulate_on_the_undeformed_configuration"].GetBool()):
            index=0
            self.curr_coords = []
            for node in self.optimization_model_part.Nodes:     
                self.curr_coords.append(node.X)            
                node.X = self.init_coords[3*index+0]
                node.X0 = node.X
                self.curr_coords.append(node.Y)
                node.Y = self.init_coords[3*index+1]
                node.Y0 = node.Y
                self.curr_coords.append(node.Z)
                node.Z = self.init_coords[3*index+2]
                node.Z0 = node.Z           
                index = index + 1   

        self.model_part_controller.ComputeUnitSurfaceNormals()

        self.mapper.InverseMap(KSO.DF1DX, KSO.DF1DX_MAPPED)

        for constraint in self.constraints:
            con_id = constraint["identifier"].GetString()
            gradient_variable = self.constraint_gradient_variables[con_id]["gradient"]
            mapped_gradient_variable = self.constraint_gradient_variables[con_id]["mapped_gradient"]
            self.mapper.InverseMap(gradient_variable, mapped_gradient_variable)

        self.__computeControlPointUpdate()

        self.mapper.Map(KSO.CONTROL_POINT_UPDATE, KSO.SHAPE_UPDATE)    

        if(self.mapper.implicit_settings["formulate_on_the_undeformed_configuration"].GetBool()):
            index=0
            for node in self.optimization_model_part.Nodes: 
                node_curr_coords = self.curr_coords[3*index:3*index+3]
                node.X = node_curr_coords[0]
                node.X0 = node.X   
                node.Y = node_curr_coords[1]
                node.Y0 = node.Y   
                node.Z = node_curr_coords[2]
                node.Z0 = node.Z                         
                index = index + 1        

    # --------------------------------------------------------------------------
    def __computeControlPointUpdate(self):
        """adapted from https://msulaiman.org/onewebmedia/GradProj_2.pdf"""
        g_a, g_a_variables = self.__getActiveConstraints()

        KM.Logger.PrintInfo("ShapeOpt", "Assemble vector of objective gradient.")
        nabla_f = KM.Vector()
        s = KM.Vector()
        self.optimization_utilities.AssembleVector(self.optimization_model_part, nabla_f, KSO.DF1DX_MAPPED)

        if len(g_a) == 0:
            KM.Logger.PrintInfo("ShapeOpt", "No constraints active, use negative objective gradient as search direction.")
            s = nabla_f * (-1.0)
            s *= self.step_size / s.norm_inf()
            self.optimization_utilities.AssignVectorToVariable(self.optimization_model_part, s, KSO.SEARCH_DIRECTION)
            self.optimization_utilities.AssignVectorToVariable(self.optimization_model_part, [0.0]*len(s), KSO.CORRECTION)
            self.optimization_utilities.AssignVectorToVariable(self.optimization_model_part, s, KSO.CONTROL_POINT_UPDATE)
            return


        KM.Logger.PrintInfo("ShapeOpt", "Assemble matrix of constraint gradient.")
        N = KM.Matrix()
        self.optimization_utilities.AssembleMatrix(self.optimization_model_part, N, g_a_variables)  # TODO check if gradients are 0.0! - in cpp

        settings = KM.Parameters('{ "solver_type" : "LinearSolversApplication.dense_col_piv_householder_qr" }')
        solver = dense_linear_solver_factory.ConstructSolver(settings)

        KM.Logger.PrintInfo("ShapeOpt", "Calculate projected search direction and correction.")
        c = KM.Vector()
        self.optimization_utilities.CalculateProjectedSearchDirectionAndCorrection(
            nabla_f,
            N,
            g_a,
            solver,
            s,
            c)

        if c.norm_inf() != 0.0:
            if c.norm_inf() <= self.max_correction_share * self.step_size:
                delta = self.step_size - c.norm_inf()
                s *= delta/s.norm_inf()
            else:
                KM.Logger.PrintWarning("ShapeOpt", f"Correction is scaled down from {c.norm_inf()} to {self.max_correction_share * self.step_size}.")
                c *= self.max_correction_share * self.step_size / c.norm_inf()
                s *= (1.0 - self.max_correction_share) * self.step_size / s.norm_inf()
        else:
            s *= self.step_size / s.norm_inf()

        self.optimization_utilities.AssignVectorToVariable(self.optimization_model_part, s, KSO.SEARCH_DIRECTION)
        self.optimization_utilities.AssignVectorToVariable(self.optimization_model_part, c, KSO.CORRECTION)
        self.optimization_utilities.AssignVectorToVariable(self.optimization_model_part, s+c, KSO.CONTROL_POINT_UPDATE)

    # --------------------------------------------------------------------------
    def __getActiveConstraints(self):
        active_constraint_values = []
        active_constraint_variables = []

        for constraint in self.constraints:
            if self.__isConstraintActive(constraint):
                identifier = constraint["identifier"].GetString()
                constraint_value = self.communicator.getStandardizedValue(identifier)
                active_constraint_values.append(constraint_value)
                active_constraint_variables.append(
                    self.constraint_gradient_variables[identifier]["mapped_gradient"])

        return active_constraint_values, active_constraint_variables

    # --------------------------------------------------------------------------
    def __isConstraintActive(self, constraint):
        identifier = constraint["identifier"].GetString()
        constraint_value = self.communicator.getStandardizedValue(identifier)
        if constraint["type"].GetString() == "=":
            return True
        elif constraint_value >= 0:
            return True
        else:
            return False

    # --------------------------------------------------------------------------
    def __logCurrentOptimizationStep(self):
        additional_values_to_log = {}
        additional_values_to_log["step_size"] = self.step_size
        additional_values_to_log["inf_norm_s"] = self.optimization_utilities.ComputeMaxNormOfNodalVariable(self.design_surface, KSO.SEARCH_DIRECTION)
        additional_values_to_log["inf_norm_c"] = self.optimization_utilities.ComputeMaxNormOfNodalVariable(self.design_surface, KSO.CORRECTION)
        self.data_logger.LogCurrentValues(self.optimization_iteration, additional_values_to_log)
        self.data_logger.LogCurrentDesign(self.optimization_iteration)

    # --------------------------------------------------------------------------
    def __isAlgorithmConverged(self):

        if self.optimization_iteration > 1 :

            # Check if maximum iterations were reached
            if self.optimization_iteration == self.max_iterations:
                KM.Logger.Print("")
                KM.Logger.PrintInfo("ShapeOpt", "Maximal iterations of optimization problem reached!")
                return True

            # Check for relative tolerance
            relative_change_of_objective_value = self.data_logger.GetValues("rel_change_objective")[self.optimization_iteration]
            if abs(relative_change_of_objective_value) < self.relative_tolerance:
                KM.Logger.Print("")
                KM.Logger.PrintInfo("ShapeOpt", "Optimization problem converged within a relative objective tolerance of ",self.relative_tolerance,"%.")
                return True

    # --------------------------------------------------------------------------
    def __determineAbsoluteChanges(self):
        self.optimization_utilities.AddFirstVariableToSecondVariable(self.design_surface, KSO.CONTROL_POINT_UPDATE, KSO.CONTROL_POINT_CHANGE)
        self.optimization_utilities.AddFirstVariableToSecondVariable(self.design_surface, KSO.SHAPE_UPDATE, KSO.SHAPE_CHANGE)

# ==============================================================================
