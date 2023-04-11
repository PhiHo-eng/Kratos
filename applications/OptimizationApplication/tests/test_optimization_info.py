
import KratosMultiphysics as Kratos

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.OptimizationApplication.utilities.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.responses.mass_response_function import MassResponseFunction

class TestOptimizationData(kratos_unittest.TestCase):
    def test_AdvanceStep(self):
        optimization_info = OptimizationInfo()
        problem_data = optimization_info.GetProblemDataContainer()
        problem_data["test"] = OptimizationInfo.OptimizationData(3)
        optimization_data = problem_data["test"]

        optimization_data["step"] = 1
        self.assertEqual(optimization_data["step"], 1)

        optimization_info.AdvanceStep()
        optimization_data["step"] = 2
        self.assertEqual(optimization_data["step"], 2)
        self.assertEqual(optimization_data["step", 1], 1)

        optimization_info.AdvanceStep()
        optimization_data["step"] = 3
        self.assertEqual(optimization_data["step"], 3)
        self.assertEqual(optimization_data["step", 1], 2)
        self.assertEqual(optimization_data["step", 2], 1)

        optimization_info.AdvanceStep()
        optimization_data["step"] = 4
        self.assertEqual(optimization_data["step"], 4)
        self.assertEqual(optimization_data["step", 1], 3)
        self.assertEqual(optimization_data["step", 2], 2)
        with self.assertRaises(RuntimeError):
            _ = optimization_data["step", 3]

    def test_AdvanceStepSet(self):
        optimization_info = OptimizationInfo()
        problem_data = optimization_info.GetProblemDataContainer()
        problem_data["test"] = OptimizationInfo.OptimizationData(3)
        optimization_data = problem_data["test"]
        optimization_data["test_1/test_sub_1"] = OptimizationInfo.OptimizationData(2)

        self.assertEqual(optimization_data["test_1"].GetBufferSize(), 3)
        self.assertEqual(optimization_data["test_1/test_sub_1"].GetBufferSize(), 2)

        optimization_data["test_1/int"] = 1
        optimization_data["test_1/test_sub_1/int"] = 2

        optimization_info.AdvanceStep()
        optimization_data["test_1/int"] = 3
        optimization_data["test_1/test_sub_1/int"] = 4

        optimization_info.AdvanceStep()
        optimization_data["test_1/int"] = 5
        optimization_data["test_1/test_sub_1/int"] = 6

        optimization_info.AdvanceStep()
        optimization_data["test_1/int"] = 7
        optimization_data["test_1/test_sub_1/int"] = 8

        optimization_info.AdvanceStep()
        optimization_data["test_1/int"] = 9
        optimization_data["test_1/test_sub_1/int"] = 10

        optimization_info.AdvanceStep()
        optimization_data["test_1/int"] = 11
        optimization_data["test_1/test_sub_1/int"] = 12

        self.assertEqual(optimization_data["test_1/int", 0], 11)
        self.assertEqual(optimization_data["test_1/int", 1], 9)
        self.assertEqual(optimization_data["test_1/int", 2], 7)
        with self.assertRaises(RuntimeError):
            optimization_data["test_1/int", 3]

        self.assertEqual(optimization_data["test_1/test_sub_1/int", 0], 12)
        self.assertEqual(optimization_data["test_1/test_sub_1/int", 1], 10)
        with self.assertRaises(RuntimeError):
            _ = optimization_data["test_1/test_sub_1/int", 3]

    def test_HasValue(self):
        optimization_info = OptimizationInfo()
        problem_data = optimization_info.GetProblemDataContainer()
        problem_data["test"] = OptimizationInfo.OptimizationData(2)
        optimization_data = problem_data["test"]
        optimization_data["data/sub_1"] = 1
        optimization_data["data/sub_2", 1] = 2
        optimization_data.SetValue("data/sub_3/sub_sub1", 3, 1)

        self.assertTrue(optimization_data.HasValue("data/sub_2", 1))
        self.assertTrue(optimization_data["data"].HasValue("sub_3/sub_sub1", 1))
        self.assertTrue(optimization_data.HasValue("data/sub_2", 1))
        self.assertTrue(optimization_data.HasValue("data/sub_3", 0))
        self.assertTrue(optimization_data.HasValue("data/sub_3", 1))
        self.assertFalse(optimization_data.HasValue("data/sub_2"))
        self.assertFalse(optimization_data.HasValue("data_2/sub_2"))

    def test_RemoveValue(self):
        optimization_info = OptimizationInfo()
        problem_data = optimization_info.GetProblemDataContainer()
        problem_data["test"] = OptimizationInfo.OptimizationData(3)
        optimization_data = problem_data["test"]
        optimization_data["data/sub_1"] = 1
        optimization_data["data/sub_1", 1] = 2
        optimization_data["data/sub_1", 2] = 3

        self.assertTrue(optimization_data.HasValue("data/sub_1", 0))
        self.assertTrue(optimization_data.HasValue("data/sub_1", 1))
        self.assertTrue(optimization_data.HasValue("data/sub_1", 2))

        del optimization_data["data/sub_1", 0]
        self.assertFalse(optimization_data.HasValue("data/sub_1", 0))
        self.assertTrue(optimization_data.HasValue("data/sub_1", 1))
        self.assertTrue(optimization_data.HasValue("data/sub_1", 2))

        del optimization_data["data/sub_1", 1]
        self.assertFalse(optimization_data.HasValue("data/sub_1", 0))
        self.assertFalse(optimization_data.HasValue("data/sub_1", 1))
        self.assertTrue(optimization_data.HasValue("data/sub_1", 2))

        del optimization_data["data/sub_1", 2]
        self.assertFalse(optimization_data.HasValue("data/sub_1", 0))
        self.assertFalse(optimization_data.HasValue("data/sub_1", 1))
        self.assertFalse(optimization_data.HasValue("data/sub_1", 2))
        self.assertTrue(optimization_data.HasValue("data"))

    def test_SetValue(self):
        optimization_info = OptimizationInfo()
        problem_data = optimization_info.GetProblemDataContainer()
        problem_data["test"] = OptimizationInfo.OptimizationData(3)
        optimization_data = problem_data["test"]
        optimization_data["data/sub_1", 1] = 1
        with self.assertRaises(RuntimeError):
            _ = optimization_data["data/sub_1/int"]

        with self.assertRaises(RuntimeError):
            _ = optimization_data["data"] = 2

    def test_GetMap(self):
        optimization_info = OptimizationInfo()
        problem_data = optimization_info.GetProblemDataContainer()
        problem_data["test"] = OptimizationInfo.OptimizationData(3)
        optimization_data = problem_data["test"]
        optimization_data["data/test/int"] = 1
        optimization_data["data/test/float"] = 3.0
        optimization_data["data/test/float", 1] = 4.0
        optimization_data["data/test/str", 1] = "test"
        optimization_data["data/str", 2] = "old"
        optimization_data["data/int"] = 10

        result = optimization_data.GetMap()
        for k, v in result.items():
            self.assertEqual(optimization_data[k], v)
        self.assertEqual(
            result,
            {
                "data": optimization_data["data"],
                "data/int": 10,
                "data/test": optimization_data["data/test"],
                "data/test/int": 1,
                "data/test/float": 3.0
            })

        result = optimization_data.GetMap(1)
        for k, v in result.items():
            self.assertEqual(optimization_data[k, 1], v)
        self.assertEqual(
            result,
            {
                "data": optimization_data["data"],
                "data/test": optimization_data["data/test"],
                "data/test/str": "test",
                "data/test/float": 4.0
            })

        result = optimization_data.GetMap(2)
        for k, v in result.items():
            self.assertEqual(optimization_data[k, 2], v)
        self.assertEqual(
            result,
            {
                "data": optimization_data["data"],
                "data/str": "old",
                "data/test": optimization_data["data/test"]
            })

    def test_Buffers(self):
        optimization_info = OptimizationInfo()
        problem_data = optimization_info.GetProblemDataContainer()
        problem_data["test"] = OptimizationInfo.OptimizationData(3)
        optimization_data = problem_data["test"]
        optimization_data["data/sub_1/sub_sub"] = OptimizationInfo.OptimizationData(4)
        optimization_data["data/sub_2"] = OptimizationInfo.OptimizationData(5)

        self.assertEqual(optimization_data.GetBufferSize(), 3)
        self.assertEqual(optimization_data["data"].GetBufferSize(), 3)
        self.assertEqual(optimization_data["data/sub_1"].GetBufferSize(), 3)
        self.assertEqual(optimization_data["data/sub_1/sub_sub"].GetBufferSize(), 4)
        self.assertEqual(optimization_data["data/sub_2"].GetBufferSize(), 5)

    def test_GetParent(self):
        optimization_info = OptimizationInfo()
        problem_data = optimization_info.GetProblemDataContainer()
        problem_data["test"] = OptimizationInfo.OptimizationData(3)
        optimization_data = problem_data["test"]
        optimization_data["data/sub_1/sub_sub"] = OptimizationInfo.OptimizationData(4)
        optimization_data["data/sub_2"] = OptimizationInfo.OptimizationData(5)

        self.assertEqual(optimization_data["data/sub_1/sub_sub"].GetParent(), optimization_data["data/sub_1"])
        self.assertEqual(optimization_data["data/sub_1"].GetParent(), optimization_data["data"])
        self.assertEqual(optimization_data["data/sub_2"].GetParent(), optimization_data["data"])

    def test_GetRoot(self):
        optimization_info = OptimizationInfo()
        problem_data = optimization_info.GetProblemDataContainer()
        problem_data["test"] = OptimizationInfo.OptimizationData(3)
        optimization_data = problem_data["test"]
        optimization_data["data/sub_1/sub_sub"] = OptimizationInfo.OptimizationData(4)
        optimization_data["data/sub_2"] = OptimizationInfo.OptimizationData(5)

        self.assertEqual(optimization_data["data/sub_1/sub_sub"].GetRoot(), problem_data.GetParent())
        self.assertEqual(optimization_data["data/sub_1"].GetRoot(), problem_data.GetParent())
        self.assertEqual(optimization_data["data/sub_2"].GetRoot(), problem_data.GetParent())
        self.assertEqual(optimization_data["data"].GetRoot(), problem_data.GetParent())
        self.assertEqual(optimization_data.GetRoot(), problem_data.GetParent())

class TestOptimizationInfo(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.optimization_info = OptimizationInfo()

        cls.response_function = MassResponseFunction(cls.model, Kratos.Parameters("""{"evaluated_model_part_names": ["test"]}"""), cls.optimization_info)
        cls.optimization_info.AddComponent("mass", cls.response_function)

    def test_GetComponentProblemDataContainer(self):
        mass_problem_data = self.optimization_info.GetComponentProblemDataContainer("mass")
        self.assertEqual(mass_problem_data, self.optimization_info.GetComponentProblemDataContainer(self.response_function))

    def test_AdvanceStep(self):
        mass_problem_data = self.optimization_info.GetComponentProblemDataContainer("mass")
        mass_problem_data["mass_1"] = OptimizationInfo.OptimizationData(3)

        initial_step = self.optimization_info.GetStep()

        mass_problem_data["mass_1/int"] = 1
        self.assertEqual(self.optimization_info.GetStep(), initial_step)

        self.optimization_info.AdvanceStep()
        mass_problem_data["mass_1/int"] = 2
        self.assertEqual(mass_problem_data["mass_1/int", 1], 1)
        self.assertEqual(self.optimization_info.GetStep(), initial_step + 1)
        self.assertEqual(self.optimization_info.GetComponent("mass"), self.response_function)

        self.optimization_info.AdvanceStep()
        mass_problem_data["mass_1/int"] = 3
        self.assertEqual(mass_problem_data["mass_1/int", 1], 2)
        self.assertEqual(mass_problem_data["mass_1/int", 2], 1)
        self.assertEqual(self.optimization_info.GetStep(), initial_step + 2)
        self.assertEqual(self.optimization_info.GetComponent("mass"), self.response_function)

        self.optimization_info.AdvanceStep()
        mass_problem_data["mass_1/int"] = 4
        self.assertEqual(mass_problem_data["mass_1/int", 1], 3)
        self.assertEqual(mass_problem_data["mass_1/int", 2], 2)
        self.assertEqual(self.optimization_info.GetStep(), initial_step + 3)
        self.assertEqual(self.optimization_info.GetComponent("mass"), self.response_function)

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()