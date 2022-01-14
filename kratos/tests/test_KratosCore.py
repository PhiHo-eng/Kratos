# import Kratos
import KratosMultiphysics

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests or test_classes to create the suites
import test_bounding_box
import test_calculate_distance_to_skin
import test_embedded_skin_mapping
import test_model_part
import test_model_part_io
import test_kratos_parameters
import test_materials_input
import test_geometries
import test_linear_solvers
import test_eigen_solvers
import test_condition_number
import test_point
import test_processes
import test_properties
import test_importing
import test_connectivity_preserve_modeler
import test_model
import test_redistance
import test_levelset_convection
import test_variable_utils
import test_reorder
import test_exact_integration
import test_gid_io
import test_vtk_output_process
import test_vector_interface
import test_matrix_interface
import test_restart
import test_gid_io_gauss_points
import test_mortar_mapper
import test_normal_utils
import test_skin_detection_process
import test_sparse_multiplication
import test_variable_component
import test_variable_redistribution
import test_object_printing
import test_array_1d_interface
import test_flags
import test_time_discretization
import test_serializer
import test_file_logger
import test_dofs
import test_time_averaging
import test_scipy_conversion_tools
import test_numpy_export_dense_matrix
import test_linear_constraints
import test_specifications_utilities
import test_cad_json_input
import test_compare_elements_conditions
import test_matrix_market_interface
import test_factories
import test_coordinate_transformation_utils
import test_sensitivity_utilities
import test_file_name_data_collector
import test_function_parser_utility
import test_integration_points
import test_mls_shape_functions_utility
import test_model_part_combination_utilities
import test_force_and_torque_utils
import test_print_info_in_file

def AssembleTestSuites():
    ''' Populates the test suites to run.

    Populates the test suites to run. At least, it should pupulate the suites:
    "small", "nighlty" and "all"

    Return
    ------

    suites: A dictionary of suites
        The set of suites with its test_cases added.
    '''

    suites = KratosUnittest.KratosSuites

    # Create a test suite with the selected tests (Small tests):
    smallSuite = suites['small']
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_bounding_box.TestBoundingBox]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_print_info_in_file.TestPrintInfoInFile]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_calculate_distance_to_skin.TestCalculateDistanceToSkin]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_embedded_skin_mapping.TestEmbeddedSkinMapping]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_model_part.TestModelPart]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_model_part_io.TestModelPartIO]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_model_part_io.TestModelPartIOMPI]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_materials_input.TestMaterialsInput]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_geometries.TestGeometry]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_kratos_parameters.TestParameters]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_linear_solvers.TestLinearSolvers]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_eigen_solvers.TestEigenSolvers]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_condition_number.TestConditionNumber]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_point.TestPoint]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_processes.TestProcesses]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_properties.TestProperties]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_importing.TestImporting]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_connectivity_preserve_modeler.TestConnectivityPreserveModeler]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_model.TestModel]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_redistance.TestRedistance]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_levelset_convection.TestLevelSetConvection]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_variable_utils.TestVariableUtils]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_reorder.TestReorder]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_exact_integration.TestExactIntegration]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_gid_io.TestGidIO]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_vtk_output_process.TestVtkOutputProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_vector_interface.TestVectorInterface]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_matrix_interface.TestMatrixInterface]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_restart.TestRestart]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_gid_io_gauss_points.TestGiDIOGaussPoints]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_skin_detection_process.TestSkinDetectionProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_normal_utils.TestNormalUtilsCoarseSphere]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_normal_utils.TestNormalUtilsQuadSphere]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_normal_utils.TestNormalUtils2DSymmetricalSquare]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_sparse_multiplication.TestSparseMatrixSum]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_sparse_multiplication.TestSparseMatrixTranspose]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_sparse_multiplication.TestSparseMatrixMultiplication]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_variable_component.TestVariableComponent]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_variable_redistribution.TestVariableRedistributionUtility]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_object_printing.TestObjectPrinting]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_array_1d_interface.TestArray1DInterface]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_flags.TestFlags]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_time_discretization.TestTimeDiscretization]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_serializer.TestSerializer]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_file_logger.TestFileLogger]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_dofs.TestDofs]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_time_averaging.TimeAveragingProcessTests]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_scipy_conversion_tools.TestScipyConversionTools]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_numpy_export_dense_matrix.TestNumpyExportDenseMatrix]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_linear_constraints.TestLinearMultipointConstraints]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_linear_constraints.TestLinearConstraints]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_specifications_utilities.TestSpecificationsUtilities]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_cad_json_input.TestCadJsonInput]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_compare_elements_conditions.TestCompareElementsAndConditionsUtility]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_matrix_market_interface.TestMatrixMarketInterface]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_factories.TestFactories]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_coordinate_transformation_utils.TestCoordinateTransformationUtilitiesCoarseSphere]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_coordinate_transformation_utils.TestCoordinateTransformationUtilities2DSymmetricalSquare]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_sensitivity_utilities.TestSensitivityUtilitiesTwoDimSymmetricalSquare]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_file_name_data_collector.TestFileNameDataCollector]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_function_parser_utility.TestGenericFunctionUtility]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_integration_points.TestIntegrationPoints]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_mls_shape_functions_utility.TestMLSShapeFunctionsUtility]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_model_part_combination_utilities.TestModelPartCombinationUtilities]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_force_and_torque_utils.TestForceAndTorqueUtils]))

    # Create a test suite with the selected tests plus all small tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_calculate_distance_to_skin.TestCalculateDistanceToSkin]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_mortar_mapper.TestMortarMapperCore]))

    # Create a test suite that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(nightSuite) # already contains the smallSuite

    return suites

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.runTests(AssembleTestSuites())
