//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// System includes


// External includes


// Project includes
#include "containers/model.h"
#include "geometries/quadrilateral_2d_4.h"
#include "processes/structured_mesh_generator_process.h"
#include "testing/testing.h"

// Application includes
#include "custom_utilities/fluid_mesh_utilities.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(AssignNeighbourElementsToConditions, FluidDynamicsApplicationFastSuite)
{
    // Set the test model part
    Model test_model;
    auto& r_test_model_part = test_model.CreateModelPart("TestModelPart");

    auto p_point_1 = Kratos::make_intrusive<Node>(1, -1.0, -1.0,  0.0);
    auto p_point_2 = Kratos::make_intrusive<Node>(2, -1.0,  1.0,  0.0);
    auto p_point_3 = Kratos::make_intrusive<Node>(3,  1.0,  1.0,  0.0);
    auto p_point_4 = Kratos::make_intrusive<Node>(4,  1.0, -1.0,  0.0);
    Quadrilateral2D4<Node> geometry(p_point_1, p_point_2, p_point_3, p_point_4);

    Parameters mesher_parameters(R"(
    {
        "number_of_divisions": 2,
        "element_name": "Element2D3N",
        "condition_name" : "Condition",
        "create_skin_sub_model_part": true
    })");
    StructuredMeshGeneratorProcess(geometry, r_test_model_part, mesher_parameters).Execute();

    // Call the neighbours assignment function
    const bool check_repeated_conditions = true;
    FluidMeshUtilities::AssignNeighbourElementsToConditions(r_test_model_part, check_repeated_conditions);

    // Collect obtained results in a matrix
    std::size_t aux_i = 0;
    Matrix results(r_test_model_part.NumberOfConditions(), 2);
    for (const auto& r_cond : r_test_model_part.Conditions()) {
        results(aux_i, 0) = r_cond.Id();
        results(aux_i, 1) = r_cond.GetValue(NEIGHBOUR_ELEMENTS)[0].Id();
        ++aux_i;
    }

    // Check neihgbours assignment
    Matrix expected_results(r_test_model_part.NumberOfConditions(), 2);
    expected_results(0,0) = 1; expected_results(0,1) = 7;
    expected_results(1,0) = 2; expected_results(1,1) = 6;
    expected_results(2,0) = 3; expected_results(2,1) = 8;
    expected_results(3,0) = 4; expected_results(3,1) = 2;
    expected_results(4,0) = 5; expected_results(4,1) = 3;
    expected_results(5,0) = 6; expected_results(5,1) = 3;
    expected_results(6,0) = 7; expected_results(6,1) = 1;
    expected_results(7,0) = 8; expected_results(7,1) = 6;
    KRATOS_CHECK_MATRIX_EQUAL(results, expected_results);
}

}