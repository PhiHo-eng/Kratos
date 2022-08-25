import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as SMA
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestSetMovingLoadProcess(KratosUnittest.TestCase):

    def _TestSetMovingLoad(self):
        """
        Tests a moving load on 1 condition element, where the nodes of the element are sorted in the direction of the
        moving load
        Returns
        -------

        """

        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)

        # create nodes
        second_coord = [1, 0, 0.0]
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,second_coord[0],second_coord[1],0.0)

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)

        # create condition
        cond = mp.CreateNewCondition("MovingLoadCondition2D2N", 1, [1,2], mp.GetProperties()[1])

        parameters = KratosMultiphysics.Parameters("""
                {
                    "help"            : "This process applies a moving load condition belonging to a modelpart. The load moves over line elements.",
                    "model_part_name" : "please_specify_model_part_name",
                    "variable_name"   : "POINT_LOAD",
                    "load"            : [0.0, -2.0, 0.0],
                    "direction"       : [1,1,1],
                    "velocity"        : 1
                }
                """
                                                         )
        mp.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0)
        mp.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 0.25)

        process = SMA.SetMovingLoadProcess(mp, parameters)

        # initialise and set load
        process.ExecuteInitialize()
        process.ExecuteInitializeSolutionStep()

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0, 0)
        rhs = KratosMultiphysics.Vector(0)

        # set load on node
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

        self.assertAlmostEqual(rhs[0], 0)
        self.assertAlmostEqual(rhs[1], -2)
        self.assertAlmostEqual(rhs[2], 0)
        self.assertAlmostEqual(rhs[3], 0)

        # move load
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        # check if interpolation is done correctly
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

        self.assertAlmostEqual(rhs[0], 0)
        self.assertAlmostEqual(rhs[1], -1.5)
        self.assertAlmostEqual(rhs[2], 0)
        self.assertAlmostEqual(rhs[3], -0.5)

    def _TestSetMovingLoadReverseGeom(self):
        """
        Tests a moving load on 1 condition element, where the nodes of the element are reversed compared to the
        direction of the moving load
        Returns
        -------

        """

        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)

        # create nodes
        second_coord = [1, 0, 0.0]
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, second_coord[0], second_coord[1], 0.0)

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z, mp)

        # create condition
        cond = mp.CreateNewCondition("MovingLoadCondition2D2N", 1, [2, 1], mp.GetProperties()[1])

        parameters = KratosMultiphysics.Parameters("""
                   {
                       "help"            : "This process applies a moving load condition belonging to a modelpart. The load moves over line elements.",
                       "model_part_name" : "please_specify_model_part_name",
                       "variable_name"   : "POINT_LOAD",
                       "load"            : [0.0, -2.0, 0.0],
                       "direction"       : [1,1,1],
                       "velocity"        : 1
                   }
                   """
                                                   )
        mp.ProcessInfo.SetValue(KratosMultiphysics.TIME,
                                0)
        mp.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME,
                                0.25)

        process = SMA.SetMovingLoadProcess(mp, parameters)

        # initialize and set load
        process.ExecuteInitialize()
        process.ExecuteInitializeSolutionStep()

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0, 0)
        rhs = KratosMultiphysics.Vector(0)

        # set load on node
        # cond.SetValue(SMA.MOVING_LOAD_LOCAL_DISTANCE, 0)
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

        self.assertAlmostEqual(rhs[0], 0)
        self.assertAlmostEqual(rhs[1], 0)
        self.assertAlmostEqual(rhs[2], 0)
        self.assertAlmostEqual(rhs[3], -2)

        # move load
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        # check if interpolation is done correctly
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

        self.assertAlmostEqual(rhs[0], 0)
        self.assertAlmostEqual(rhs[1], -0.5)
        self.assertAlmostEqual(rhs[2], 0)
        self.assertAlmostEqual(rhs[3], -1.5)

    def _TestSetMovingLoadMultipleConditions(self):
        """
        Tests a moving load on 2 condition elements, where the order of the elements is sorted in the direction of the
        moving load
        Returns
        -------

        """

        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)

        #create nodes
        second_coord = [1.0, 0.0, 0.0]
        third_coord = [2.0, 0.0, 0.0]
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, second_coord[0],second_coord[1],second_coord[2])
        mp.CreateNewNode(3, third_coord[0], third_coord[1], third_coord[2])


        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)

        # create condition
        conditions=[]
        conditions.append(mp.CreateNewCondition("MovingLoadCondition2D2N", 1, [1, 2], mp.GetProperties()[1]))
        conditions.append(mp.CreateNewCondition("MovingLoadCondition2D2N", 2, [2, 3], mp.GetProperties()[1]))

        parameters = KratosMultiphysics.Parameters("""
                {
                    "help"            : "This process applies a moving load condition belonging to a modelpart. The load moves over line elements.",
                    "model_part_name" : "please_specify_model_part_name",
                    "variable_name"   : "POINT_LOAD",
                    "load"            : [0.0, -2.0, 0.0],
                    "direction"       : [1,1,1],
                    "velocity"        : 1
                }
                """
                                                         )
        mp.ProcessInfo.SetValue(KratosMultiphysics.TIME,
                                                  0)
        mp.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME,
                                                  0.5)
        process = SMA.SetMovingLoadProcess(mp,parameters)


        # initialize and set load
        process.ExecuteInitialize()
        process.ExecuteInitializeSolutionStep()

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0,0)
        rhs = KratosMultiphysics.Vector(0)

        # set load on node

        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.assertAlmostEqual(all_rhs[0][0], 0)
        self.assertAlmostEqual(all_rhs[0][1], -2)
        self.assertAlmostEqual(all_rhs[0][2], 0)
        self.assertAlmostEqual(all_rhs[0][3], 0)

        self.assertAlmostEqual(all_rhs[1][0], 0)
        self.assertAlmostEqual(all_rhs[1][1], 0)
        self.assertAlmostEqual(all_rhs[1][2], 0)
        self.assertAlmostEqual(all_rhs[1][3], 0)

        # move load within first element
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        # check if interpolation is done correctly
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.assertAlmostEqual(all_rhs[0][0], 0)
        self.assertAlmostEqual(all_rhs[0][1], -1)
        self.assertAlmostEqual(all_rhs[0][2], 0)
        self.assertAlmostEqual(all_rhs[0][3], -1)

        self.assertAlmostEqual(all_rhs[1][0], 0)
        self.assertAlmostEqual(all_rhs[1][1], 0)
        self.assertAlmostEqual(all_rhs[1][2], 0)
        self.assertAlmostEqual(all_rhs[1][3], 0)

        # move load to element connection element
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.assertAlmostEqual(all_rhs[0][0], 0)
        self.assertAlmostEqual(all_rhs[0][1], 0)
        self.assertAlmostEqual(all_rhs[0][2], 0)
        self.assertAlmostEqual(all_rhs[0][3], -2)

        self.assertAlmostEqual(all_rhs[1][0], 0)
        self.assertAlmostEqual(all_rhs[1][1], 0)
        self.assertAlmostEqual(all_rhs[1][2], 0)
        self.assertAlmostEqual(all_rhs[1][3], 0)

        # move load to next element
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.assertAlmostEqual(all_rhs[0][0], 0)
        self.assertAlmostEqual(all_rhs[0][1], 0)
        self.assertAlmostEqual(all_rhs[0][2], 0)
        self.assertAlmostEqual(all_rhs[0][3], 0)

        self.assertAlmostEqual(all_rhs[1][0], 0)
        self.assertAlmostEqual(all_rhs[1][1], -1)
        self.assertAlmostEqual(all_rhs[1][2], 0)
        self.assertAlmostEqual(all_rhs[1][3], -1)


    def _TestSetMovingLoadMultipleConditionsReversed(self):
        """
        Tests a moving load on 2 condition elements, where the order of the elements is reversed compared to the moving
        direction of the load
        Returns
        -------

        """

        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)

        # create nodes
        second_coord = [1.0, 0.0, 0.0]
        third_coord = [2.0, 0.0, 0.0]
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, second_coord[0],second_coord[1],second_coord[2])
        mp.CreateNewNode(3, third_coord[0], third_coord[1], third_coord[2])

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)

        # create condition
        conditions = []
        conditions.append(mp.CreateNewCondition("MovingLoadCondition2D2N", 1, [3, 2], mp.GetProperties()[1]))
        conditions.append(mp.CreateNewCondition("MovingLoadCondition2D2N", 2, [2, 1], mp.GetProperties()[1]))

        # set parameters and process info
        parameters = KratosMultiphysics.Parameters("""
                {
                    "help"            : "This process applies a moving load condition belonging to a modelpart. The load moves over line elements.",
                    "model_part_name" : "please_specify_model_part_name",
                    "variable_name"   : "POINT_LOAD",
                    "load"            : [0.0, -2.0, 0.0],
                    "direction"       : [1,1,1],
                    "velocity"        : 1
                }
                """
                                                         )
        mp.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0)
        mp.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 0.5)
        process = SMA.SetMovingLoadProcess(mp,parameters)

        # initialize and set load
        process.ExecuteInitialize()
        process.ExecuteInitializeSolutionStep()

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0,0)
        rhs = KratosMultiphysics.Vector(0)

        # calculate load
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.assertAlmostEqual(all_rhs[0][0], 0)
        self.assertAlmostEqual(all_rhs[0][1], 0)
        self.assertAlmostEqual(all_rhs[0][2], 0)
        self.assertAlmostEqual(all_rhs[0][3], 0)

        self.assertAlmostEqual(all_rhs[1][0], 0)
        self.assertAlmostEqual(all_rhs[1][1], 0)
        self.assertAlmostEqual(all_rhs[1][2], 0)
        self.assertAlmostEqual(all_rhs[1][3], -2)

        # move load within first element
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        # check if interpolation is done correctly
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.assertAlmostEqual(all_rhs[0][0], 0)
        self.assertAlmostEqual(all_rhs[0][1], 0)
        self.assertAlmostEqual(all_rhs[0][2], 0)
        self.assertAlmostEqual(all_rhs[0][3], 0)

        self.assertAlmostEqual(all_rhs[1][0], 0)
        self.assertAlmostEqual(all_rhs[1][1], -1)
        self.assertAlmostEqual(all_rhs[1][2], 0)
        self.assertAlmostEqual(all_rhs[1][3], -1)

        # move load to element connection element
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        # calculate load
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.assertAlmostEqual(all_rhs[0][0], 0)
        self.assertAlmostEqual(all_rhs[0][1], 0)
        self.assertAlmostEqual(all_rhs[0][2], 0)
        self.assertAlmostEqual(all_rhs[0][3], 0)

        self.assertAlmostEqual(all_rhs[1][0], 0)
        self.assertAlmostEqual(all_rhs[1][1], -2)
        self.assertAlmostEqual(all_rhs[1][2], 0)
        self.assertAlmostEqual(all_rhs[1][3], 0)

        # move load to next element
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteInitializeSolutionStep()

        # calculate load
        all_rhs = []
        for cond in conditions:
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            all_rhs.append(list(rhs))

        self.assertAlmostEqual(all_rhs[0][0], 0)
        self.assertAlmostEqual(all_rhs[0][1], -1)
        self.assertAlmostEqual(all_rhs[0][2], 0)
        self.assertAlmostEqual(all_rhs[0][3], -1)

        self.assertAlmostEqual(all_rhs[1][0], 0)
        self.assertAlmostEqual(all_rhs[1][1], 0)
        self.assertAlmostEqual(all_rhs[1][2], 0)
        self.assertAlmostEqual(all_rhs[1][3], 0)

    def test_SetMovingLoad(self):
        self._TestSetMovingLoad()

    def test_SetMovingLoadReverseGeom(self):
        self._TestSetMovingLoadReverseGeom()

    def test_SetMovingLoadMultipleConditions(self):
        self._TestSetMovingLoadMultipleConditions()

    def test_SetMovingLoadMultipleConditionsReversed(self):
        self._TestSetMovingLoadMultipleConditionsReversed()



if __name__ == '__main__':
    KratosUnittest.main()
