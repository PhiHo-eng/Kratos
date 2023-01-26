# We import the libraries
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestCheckSameModelPartUsingSkinDistanceProcess(KratosUnittest.TestCase):
    def tearDown(self):
        pass

    def setUp(self):
        # We create the model part
        self.current_model = KratosMultiphysics.Model()
        self.model_part_1 = self.current_model.CreateModelPart("model_part_1")
        self.model_part_1.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 3)
        self.model_part_1.AddNodalSolutionStepVariable(KratosMultiphysics.BULK_MODULUS)
        self.model_part_1.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_VAUX)
        self.model_part_1.AddNodalSolutionStepVariable(KratosMultiphysics.EXTERNAL_FORCES_VECTOR)
        self.model_part_1.AddNodalSolutionStepVariable(KratosMultiphysics.LOCAL_AXES_MATRIX)
        self.model_part_1.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)
        self.model_part_2 = self.current_model.CreateModelPart("model_part_2")
        self.model_part_2.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 3)
        self.model_part_2.AddNodalSolutionStepVariable(KratosMultiphysics.BULK_MODULUS)
        self.model_part_2.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_VAUX)
        self.model_part_2.AddNodalSolutionStepVariable(KratosMultiphysics.EXTERNAL_FORCES_VECTOR)
        self.model_part_2.AddNodalSolutionStepVariable(KratosMultiphysics.LOCAL_AXES_MATRIX)
        self.model_part_2.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)

        input_mdpa = GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/coarse_sphere_with_conditions")
        model_part_io = KratosMultiphysics.ModelPartIO(input_mdpa)
        model_part_io.ReadModelPart(self.model_part_1)
        model_part_io.ReadModelPart(self.model_part_2)

    def test_same_sphere(self):
        parameters = KratosMultiphysics.Parameters("""
        {
            "skin_model_part_1_name" : "",
            "skin_model_part_2_name" : ""
        }
        """)
        parameters["skin_model_part_1_name"].SetString(self.model_part_1.Name + ".Skin_Part")
        parameters["skin_model_part_2_name"].SetString(self.model_part_2.Name + ".Skin_Part")
        KratosMultiphysics.CheckSameModelPartUsingSkinDistanceProcess3D(self.current_model, parameters).Execute()

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
