import KratosMultiphysics as Kratos
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.kratos_utilities as KratosUtils

class FluidAuxiliaryUtilitiesTest(UnitTest.TestCase):

    def setUp(self):
        self.model = Kratos.Model()
        fluid_model_part = self.model.CreateModelPart("FluidModelPart")
        fluid_model_part.SetBufferSize(2)
        fluid_model_part.ProcessInfo.SetValue(Kratos.DOMAIN_SIZE, 2)
        fluid_model_part.AddNodalSolutionStepVariable(Kratos.NORMAL)
        fluid_model_part.AddNodalSolutionStepVariable(Kratos.DISTANCE)
        fluid_model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        fluid_model_part.AddNodalSolutionStepVariable(Kratos.MESH_VELOCITY)
        Kratos.ModelPartIO("Cavity/square5").ReadModelPart(fluid_model_part)

    def testCreateSlipMultiPointConstraints(self):
        # Add DOFs to model part
        fluid_model_part = self.model.GetModelPart("FluidModelPart")
        Kratos.VariableUtils().AddDof(Kratos.VELOCITY_X, fluid_model_part)
        Kratos.VariableUtils().AddDof(Kratos.VELOCITY_Y, fluid_model_part)

        # Set a slip BC node
        node_1 = fluid_model_part.GetNode(1)
        node_1.Set(Kratos.SLIP, True)
        node_1.SetSolutionStepValue(Kratos.NORMAL, [2.0, 1.0, 0.0])
        KratosFluid.FluidAuxiliaryUtilities.CreateSlipMultiPointConstraints(fluid_model_part, Kratos.VELOCITY, Kratos.SLIP)

        # Check the slip MPCs arrays
        mpc_1 = fluid_model_part.GetMasterSlaveConstraint(1)
        slave_dofs_vect = mpc_1.GetSlaveDofsVector()
        master_dofs_vect = mpc_1.GetMasterDofsVector()
        relation_matrix = Kratos.Matrix()
        constant_vector = Kratos.Vector()
        mpc_1.CalculateLocalSystem(relation_matrix, constant_vector, fluid_model_part.ProcessInfo)

        self.assertEqual(slave_dofs_vect[0].GetVariable(), Kratos.VELOCITY_X)
        self.assertEqual(master_dofs_vect[0].GetVariable(), Kratos.VELOCITY_Y)
        self.assertAlmostEqual(constant_vector[0], 0.0, 12)
        self.assertAlmostEqual(relation_matrix[0,0], -0.5, 12)

    def testCreateSlipMultiPointConstraintsWithWallVelocity(self):
        # Add DOFs to model part
        fluid_model_part = self.model.GetModelPart("FluidModelPart")
        Kratos.VariableUtils().AddDof(Kratos.VELOCITY_X, fluid_model_part)
        Kratos.VariableUtils().AddDof(Kratos.VELOCITY_Y, fluid_model_part)

        # Set a slip BC node
        node_1 = fluid_model_part.GetNode(1)
        node_1.Set(Kratos.SLIP, True)
        node_1.SetSolutionStepValue(Kratos.NORMAL, [2.0, 1.0, 0.0])
        node_1.SetSolutionStepValue(Kratos.MESH_VELOCITY, [1.0, 1.0, 0.0])
        KratosFluid.FluidAuxiliaryUtilities.CreateSlipMultiPointConstraints(fluid_model_part, Kratos.VELOCITY, Kratos.MESH_VELOCITY, Kratos.SLIP)

        # Check the slip MPCs arrays
        mpc_1 = fluid_model_part.GetMasterSlaveConstraint(1)
        slave_dofs_vect = mpc_1.GetSlaveDofsVector()
        master_dofs_vect = mpc_1.GetMasterDofsVector()
        relation_matrix = Kratos.Matrix()
        constant_vector = Kratos.Vector()
        mpc_1.CalculateLocalSystem(relation_matrix, constant_vector, fluid_model_part.ProcessInfo)

        self.assertEqual(slave_dofs_vect[0].GetVariable(), Kratos.VELOCITY_X)
        self.assertEqual(master_dofs_vect[0].GetVariable(), Kratos.VELOCITY_Y)
        self.assertAlmostEqual(constant_vector[0], 1.5, 12)
        self.assertAlmostEqual(relation_matrix[0,0], -0.5, 12)

    def testCalculateFluidVolume(self):
        fluid_model_part = self.model.GetModelPart("FluidModelPart")
        fluid_volume = KratosFluid.FluidAuxiliaryUtilities.CalculateFluidVolume(fluid_model_part)
        self.assertAlmostEqual(fluid_volume, 1.0, 12)

    def testCalculateFluidPositiveVolume(self):
        # Set fluid level set
        level_set_y = 1.0/3.0
        fluid_model_part = self.model.GetModelPart("FluidModelPart")
        for node in fluid_model_part.Nodes:
            node.SetSolutionStepValue(Kratos.DISTANCE, 0, node.Y - level_set_y)

        # Calculate the fluid positive volume
        fluid_positive_volume = KratosFluid.FluidAuxiliaryUtilities.CalculateFluidPositiveVolume(fluid_model_part)
        self.assertAlmostEqual(fluid_positive_volume, 1.0 - level_set_y, 12)

    def testCalculateFluidNegativeVolume(self):
        # Set fluid level set
        level_set_y = 2.0/3.0
        fluid_model_part = self.model.GetModelPart("FluidModelPart")
        for node in fluid_model_part.Nodes:
            node.SetSolutionStepValue(Kratos.DISTANCE, 0, node.Y - level_set_y)

        # Calculate the fluid negative volume
        fluid_negative_volume = KratosFluid.FluidAuxiliaryUtilities.CalculateFluidNegativeVolume(fluid_model_part)
        self.assertAlmostEqual(fluid_negative_volume, level_set_y, 12)

    def testCalculateFlowRatePositiveSkin(self):
        # Set fluid level set
        level_set_y = 2.0/3.0
        fluid_model_part = self.model.GetModelPart("FluidModelPart")
        for node in fluid_model_part.Nodes:
            node.SetSolutionStepValue(Kratos.VELOCITY, 0, [1.0,0.0,0.0])
            node.SetSolutionStepValue(Kratos.DISTANCE, 0, node.Y - level_set_y)

        # Call the tetrahedral mesh orientation process to calculate the normals and neighbours
        tmoc = Kratos.TetrahedralMeshOrientationCheck
        throw_errors = False
        flags = (tmoc.COMPUTE_NODAL_NORMALS).AsFalse() | (tmoc.COMPUTE_CONDITION_NORMALS).AsFalse() | tmoc.ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS
        Kratos.TetrahedralMeshOrientationCheck(fluid_model_part, throw_errors, flags).Execute()

        # Calculate the left wall flow
        ref_positive_flow_rate_left = -1.0/3.0
        left_skin_model_part = self.model.GetModelPart("FluidModelPart.NoSlip2D_left_wall")
        for cond in left_skin_model_part.Conditions:
            cond.Set(Kratos.INLET, True)
        positive_flow_rate_left = KratosFluid.FluidAuxiliaryUtilities.CalculateFlowRatePositiveSkin(left_skin_model_part, Kratos.INLET)
        self.assertAlmostEqual(positive_flow_rate_left, ref_positive_flow_rate_left, 12)

        # Calculate the right wall flow
        ref_positive_flow_rate_right = 1.0/3.0
        right_skin_model_part = self.model.GetModelPart("FluidModelPart.NoSlip2D_right_wall")
        for cond in right_skin_model_part.Conditions:
            cond.Set(Kratos.OUTLET, True)
        positive_flow_rate_right = KratosFluid.FluidAuxiliaryUtilities.CalculateFlowRatePositiveSkin(right_skin_model_part, Kratos.OUTLET)
        self.assertAlmostEqual(positive_flow_rate_right, ref_positive_flow_rate_right, 12)

        # Calculate the top wall flow (without flag)
        ref_positive_flow_rate_top = 0.0
        top_skin_model_part = self.model.GetModelPart("FluidModelPart.NoSlip2D_top_wall")
        positive_flow_rate_top = KratosFluid.FluidAuxiliaryUtilities.CalculateFlowRatePositiveSkin(top_skin_model_part)
        self.assertAlmostEqual(positive_flow_rate_top, ref_positive_flow_rate_top, 12)

    def testCalculateFlowRateNegativeSkin(self):
        # Set fluid level set
        level_set_y = 2.0/3.0
        fluid_model_part = self.model.GetModelPart("FluidModelPart")
        for node in fluid_model_part.Nodes:
            node.SetSolutionStepValue(Kratos.VELOCITY, 0, [1.0,0.0,0.0])
            node.SetSolutionStepValue(Kratos.DISTANCE, 0, node.Y - level_set_y)

        # Call the tetrahedral mesh orientation process to calculate the normals and neighbours
        tmoc = Kratos.TetrahedralMeshOrientationCheck
        throw_errors = False
        flags = (tmoc.COMPUTE_NODAL_NORMALS).AsFalse() | (tmoc.COMPUTE_CONDITION_NORMALS).AsFalse() | tmoc.ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS
        Kratos.TetrahedralMeshOrientationCheck(fluid_model_part, throw_errors, flags).Execute()

        # Calculate the left wall flow
        ref_negative_flow_rate_left = -2.0/3.0
        left_skin_model_part = self.model.GetModelPart("FluidModelPart.NoSlip2D_left_wall")
        for cond in left_skin_model_part.Conditions:
            cond.Set(Kratos.INLET, True)
        negative_flow_rate_left = KratosFluid.FluidAuxiliaryUtilities.CalculateFlowRateNegativeSkin(left_skin_model_part, Kratos.INLET)
        self.assertAlmostEqual(negative_flow_rate_left, ref_negative_flow_rate_left, 12)

        # Calculate the right wall flow
        ref_negative_flow_rate_right = 2.0/3.0
        right_skin_model_part = self.model.GetModelPart("FluidModelPart.NoSlip2D_right_wall")
        for cond in right_skin_model_part.Conditions:
            cond.Set(Kratos.OUTLET, True)
        negative_flow_rate_right = KratosFluid.FluidAuxiliaryUtilities.CalculateFlowRateNegativeSkin(right_skin_model_part, Kratos.OUTLET)
        self.assertAlmostEqual(negative_flow_rate_right, ref_negative_flow_rate_right, 12)

        # Calculate the top wall flow (without flag)
        ref_negative_flow_rate_top = 0.0
        top_skin_model_part = self.model.GetModelPart("FluidModelPart.NoSlip2D_top_wall")
        negative_flow_rate_top = KratosFluid.FluidAuxiliaryUtilities.CalculateFlowRateNegativeSkin(top_skin_model_part)
        self.assertAlmostEqual(negative_flow_rate_top, ref_negative_flow_rate_top, 12)

    def tearDown(self):
        KratosUtils.DeleteFileIfExisting("Cavity/square5.time")


if __name__ == '__main__':
    UnitTest.main()
