import KratosMultiphysics
import KratosMultiphysics.RomApplication as romapp
import KratosMultiphysics.StructuralMechanicsApplication
from KratosMultiphysics.RomApplication.empirical_cubature_method import EmpiricalCubatureMethod
from KratosMultiphysics.RomApplication import python_solvers_wrapper_rom as solver_wrapper
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
from KratosMultiphysics.RomApplication.randomized_singular_value_decomposition import RandomizedSingularValueDecomposition#####

import json
import numpy as np

class StructuralMechanicsAnalysisROM(StructuralMechanicsAnalysis):

    def __init__(self,model,project_parameters, hyper_reduction_element_selector = None, build_petrov_galerkin=False, solve_petrov_galerkin=False, solve_least_squares=False):
        self.build_petrov_galerkin = build_petrov_galerkin
        self.solve_petrov_galerkin = solve_petrov_galerkin
        self.solve_least_squares = solve_least_squares
        super().__init__(model,project_parameters)
        if hyper_reduction_element_selector != None :
            if hyper_reduction_element_selector == "EmpiricalCubature":
                self.hyper_reduction_element_selector = EmpiricalCubatureMethod()
                self.time_step_residual_matrix_container = []
            else:
                err_msg =  "The requested element selection method \"" + hyper_reduction_element_selector + "\" is not in the rom application\n"
                err_msg += "Available options are: \"EmpiricalCubature\""
                raise Exception(err_msg)
        else:
            self.hyper_reduction_element_selector = None

    #### Internal functions ####
    def _CreateSolver(self):
        """ Create the Solver (and create and import the ModelPart if it is not alread in the model) """
        ## Solver construction
        with open('RomParameters.json') as rom_parameters:
            rom_settings = KratosMultiphysics.Parameters(rom_parameters.read())
            self.project_parameters["solver_settings"].AddValue("rom_settings", rom_settings["rom_settings"])
            if self.build_petrov_galerkin:
                self.project_parameters["solver_settings"].AddBool("build_petrov_galerkin", self.build_petrov_galerkin)
            if self.solve_petrov_galerkin:
                self.project_parameters["solver_settings"].AddBool("solve_petrov_galerkin", self.solve_petrov_galerkin)
                self.project_parameters["solver_settings"].AddValue("rom_residual_settings",rom_settings["Petrov_Galerkin_basis"]["rom_settings"])
            if self.solve_least_squares:
                self.project_parameters["solver_settings"].AddBool("solve_least_squares", self.solve_least_squares)
        return solver_wrapper.CreateSolverByParameters(self.model, self.project_parameters["solver_settings"],self.project_parameters["problem_data"]["parallel_type"].GetString())

    def _GetSimulationName(self):
        return "::[ROM Simulation]:: "

    def ModifyAfterSolverInitialize(self):
        """Here is where the ROM_BASIS is imposed to each node"""
        super().ModifyAfterSolverInitialize()
        computing_model_part = self._solver.GetComputingModelPart()
        with open('RomParameters.json') as f:
            data = json.load(f)
            nodal_dofs = len(data["rom_settings"]["nodal_unknowns"])
            nodal_modes = data["nodal_modes"]
            counter = 0
            rom_dofs= self.project_parameters["solver_settings"]["rom_settings"]["number_of_rom_dofs"].GetInt()
            if self.solve_petrov_galerkin:
                nodal_residual_modes = data["Petrov_Galerkin_basis"]["nodal_modes"]####ADDEDPETROV
                nodal_residual_dofs = len(data["Petrov_Galerkin_basis"]["rom_settings"]["nodal_unknowns"])####ADDEDPETROV
                rom_residual_dofs = data["Petrov_Galerkin_basis"]["rom_settings"]["number_of_rom_dofs"]
                for node in computing_model_part.Nodes:
                    aux = KratosMultiphysics.Matrix(nodal_dofs, rom_dofs)
                    aux_residual = KratosMultiphysics.Matrix(nodal_residual_dofs,rom_residual_dofs)####ADDEDPETROV
                    for j in range(nodal_dofs):
                        Counter=str(node.Id)
                        for i in range(rom_residual_dofs):
                            if (i<rom_dofs):
                                aux[j,i] = nodal_modes[Counter][j][i]
                            aux_residual[j,i] = nodal_residual_modes[Counter][j][i]####ADDEDPETROV
                    node.SetValue(romapp.ROM_BASIS, aux ) # ROM basis
                    node.SetValue(romapp.ROM_BASIS_ASSEMBLED_RESIDUALS, aux_residual)####ADDEDPETROV
                    counter+=1
            else:
                for node in computing_model_part.Nodes:
                    aux = KratosMultiphysics.Matrix(nodal_dofs, rom_dofs)
                    for j in range(nodal_dofs):
                        Counter=str(node.Id)
                        for i in range(rom_dofs):
                            aux[j,i] = nodal_modes[Counter][j][i]
                    node.SetValue(romapp.ROM_BASIS, aux ) # ROM basis
                    counter+=1
        if self.hyper_reduction_element_selector != None:
            if self.hyper_reduction_element_selector.Name == "EmpiricalCubature":
<<<<<<< HEAD
                hyper_reduction_parameters = KratosMultiphysics.Parameters()
                hyper_reduction_parameters.AddEmptyList("nodal_unknowns")
                hyper_reduction_parameters["nodal_unknowns"].SetStringArray(self.project_parameters["solver_settings"]["rom_settings"]["nodal_unknowns"].GetStringArray())
                hyper_reduction_parameters.AddInt("number_of_rom_dofs",self.project_parameters["solver_settings"]["rom_settings"]["number_of_rom_dofs"].GetInt())
                self.ResidualUtilityObject = romapp.RomResidualsUtility(self._GetSolver().GetComputingModelPart(), hyper_reduction_parameters, self._GetSolver().get_solution_scheme())
                
=======
                self.ResidualUtilityObject = romapp.RomResidualsUtility(self._GetSolver().GetComputingModelPart(), self.project_parameters["solver_settings"]["rom_settings"], self._GetSolver()._GetScheme())
>>>>>>> origin/master

    def FinalizeSolutionStep(self):
        if self.hyper_reduction_element_selector != None:
            if self.hyper_reduction_element_selector.Name == "EmpiricalCubature":
                print('\n\n\n\nGenerating matrix of residuals')
                ResMat = self.ResidualUtilityObject.GetResiduals()
                NP_ResMat = np.array(ResMat, copy=False)
                self.time_step_residual_matrix_container.append(NP_ResMat)
        super().FinalizeSolutionStep()

    def Finalize(self):
        super().Finalize()
        if self.hyper_reduction_element_selector != None:
            if self.hyper_reduction_element_selector.Name == "EmpiricalCubature":
                OriginalNumberOfElements = self._GetSolver().GetComputingModelPart().NumberOfElements()
                ModelPartName = self._GetSolver().settings["model_import_settings"]["input_filename"].GetString()
<<<<<<< HEAD
                ###### Sebastian #######
                builder_and_solver = self._GetSolver().get_builder_and_solver()
                CalculateReactionsFlag = builder_and_solver.GetCalculateReactionsFlag()
                # Load Residual parameters to force considering some restricted elements
                if CalculateReactionsFlag:
                    try:
                        with open('RomParameters.json') as s:
                            HR_data_residuals = json.load(s)
                            restricted_residual_elements = np.array(HR_data_residuals["residual_parameters"]["rom_settings"]["restricted_elements"])-1
                        selected_restricted_residual_elements = np.random.choice(restricted_residual_elements, 5).tolist()
                        selected_restricted_residual_elements.sort()
                        self. hyper_reduction_element_selector.SetUp(self.time_step_residual_matrix_container, OriginalNumberOfElements, ModelPartName,selected_restricted_residual_elements)
                        self.hyper_reduction_element_selector.Run()
                    except:
                        self. hyper_reduction_element_selector.SetUp(self.time_step_residual_matrix_container, OriginalNumberOfElements, ModelPartName)
                        self.hyper_reduction_element_selector.Run() 
                else:
                ########################
                    self. hyper_reduction_element_selector.SetUp(self.time_step_residual_matrix_container, OriginalNumberOfElements, ModelPartName)
                    self.hyper_reduction_element_selector.Run()

    def ArrangeSnapshotMatrix(self,ResidualSnapshots):
        ### Building the Snapshot matrix ####
        for i in range (len(ResidualSnapshots)):
            if i == 0:
                SnapshotMatrix = ResidualSnapshots[i]
            else:
                SnapshotMatrix = np.c_[SnapshotMatrix,ResidualSnapshots[i]]
        return SnapshotMatrix



=======
                self. hyper_reduction_element_selector.SetUp(self.time_step_residual_matrix_container, OriginalNumberOfElements, ModelPartName)
                self.hyper_reduction_element_selector.Run()
>>>>>>> origin/master
