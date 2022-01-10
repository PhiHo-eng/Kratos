import os
import KratosMultiphysics
from KratosMultiphysics import Logger
Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
import KratosMultiphysics.DEMApplication as DEM
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.DEMApplication.DEM_analysis_stage import DEMAnalysisStage

import KratosMultiphysics.kratos_utilities as kratos_utils
import auxiliary_functions_for_tests

this_working_dir_backup = os.getcwd()

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class DEM2D_ControlModuleTestSolution(DEMAnalysisStage, KratosUnittest.TestCase):

    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM2D_control_module_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())


    def Initialize(self):
        super().Initialize()

        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM2D_control_module_tests_files")
        cm_project_parameters_file_name = os.path.join(path, "cm_parameters.json")

        with open(cm_project_parameters_file_name,'r') as parameters_file:
            self.cm_project_parameters = KratosMultiphysics.Parameters(parameters_file.read())

        #NOTE: We will transform CM utility into a process eventually
        from KratosMultiphysics.DEMApplication.multiaxial_control_module_generalized_2d_utility import MultiaxialControlModuleGeneralized2DUtility
        self.multiaxial_control_module = MultiaxialControlModuleGeneralized2DUtility(self.model, self.cm_project_parameters)
        self.multiaxial_control_module.ExecuteInitialize()

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

        self.multiaxial_control_module.ExecuteInitializeSolutionStep()

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

        self.multiaxial_control_module.ExecuteFinalizeSolutionStep()

    def PrintResultsForGid(self, time):
        super().PrintResultsForGid(time)

        self.multiaxial_control_module.PrintResults()

    def Finalize(self):
        tolerance = 1.001
        node_5_found = False
        node_6_found = False
        for node in self.rigid_face_model_part.Nodes:
            if node.Id == 5:
                node_force_x = node.GetSolutionStepValue(DEM.CONTACT_FORCES_X)
                expected_value = 316.79
                self.assertAlmostEqual(node_force_x, expected_value, delta=tolerance)
                node_5_found = True
            elif node.Id == 6:
                node_force_y = node.GetSolutionStepValue(DEM.CONTACT_FORCES_Y)
                expected_value = 150.1
                self.assertAlmostEqual(node_force_y, expected_value, delta=tolerance)
                node_6_found = True
        self.assertTrue(node_5_found)
        self.assertTrue(node_6_found)

        super().Finalize()

class DEM2D_ControlModuleTestSolutionRadial(DEMAnalysisStage, KratosUnittest.TestCase):

    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM2D_control_module_tests_files")

    def Finalize(self):
        tolerance = 1.0e-9
        first_node_found = False
        second_node_found = False
        for node in self.rigid_face_model_part.Nodes:
            if node.Id == 95:
                self.assertAlmostEqual(node.X, 0.7809379840517369, delta=tolerance)
                self.assertAlmostEqual(node.Y, 0.622777258928492, delta=tolerance)
                first_node_found = True
            elif node.Id == 110:
                self.assertAlmostEqual(node.X, -0.6227772594279198, delta=tolerance)
                self.assertAlmostEqual(node.Y, 0.7809379836521931, delta=tolerance)
                second_node_found = True
        self.assertTrue(first_node_found)
        self.assertTrue(second_node_found)

        super().Finalize()

class DEM2D_ControlModuleTestSolutionRadialMultiDofs(DEMAnalysisStage, KratosUnittest.TestCase):

    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM2D_control_module_tests_files")

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()
        if self.time <= 0.5:
            for node in self.rigid_face_model_part.Nodes:
                node.SetValue(DEM.RADIAL_NORMAL_STRESS_COMPONENT, -1.0e6 * self.time)
        else:
            for node in self.rigid_face_model_part.Nodes:
                node.SetValue(DEM.RADIAL_NORMAL_STRESS_COMPONENT, -1.0e6)

    def Finalize(self):
        tolerance = 1.0e-9
        first_node_found = False
        second_node_found = False
        for node in self.rigid_face_model_part.Nodes:
            if node.Id == 95:
                self.assertAlmostEqual(node.X, 0.636366213695056, delta=tolerance)
                self.assertAlmostEqual(node.Y, 0.5074851196038817, delta=tolerance)
                first_node_found = True
            elif node.Id == 110:
                self.assertAlmostEqual(node.X, -0.446450996247988, delta=tolerance)
                self.assertAlmostEqual(node.Y, 0.5598318428159778, delta=tolerance)
                second_node_found = True
        self.assertTrue(first_node_found)
        self.assertTrue(second_node_found)

        self.procedures.RemoveFoldersWithResults(str(self.main_path), str(self.problem_name), '')
        super().Finalize()

class TestDEM2DControlModule(KratosUnittest.TestCase):

    def setUp(self):
        pass

    def test_DEM2DControlModule(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM2D_control_module_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = KratosMultiphysics.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(DEM2D_ControlModuleTestSolution, model, parameters_file_name, auxiliary_functions_for_tests.GetHardcodedNumberOfThreads())

    def test_DEM2DControlModuleRadial(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM2D_control_module_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEMRadialCM.json")
        model = KratosMultiphysics.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(DEM2D_ControlModuleTestSolution, model, parameters_file_name, auxiliary_functions_for_tests.GetHardcodedNumberOfThreads())

    def test_DEM2DControlModuleRadialMultiDofs(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM2D_control_module_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEMRadialMultiDofsCM.json")
        model = KratosMultiphysics.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(DEM2D_ControlModuleTestSolution, model, parameters_file_name, auxiliary_functions_for_tests.GetHardcodedNumberOfThreads())

    def tearDown(self):
        file_to_remove = os.path.join("DEM2D_control_module_tests_files", "TimesPartialRelease")
        kratos_utils.DeleteFileIfExisting(GetFilePath(file_to_remove))
        os.chdir(this_working_dir_backup)

if __name__ == "__main__":
    Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
