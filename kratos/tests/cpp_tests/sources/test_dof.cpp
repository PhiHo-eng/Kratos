//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "includes/model_part.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(DofConstructorWtihoutVariableInVariablesList, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& model_part = current_model.CreateModelPart("TestModelPart");

    model_part.SetBufferSize(1);

    auto p_node = model_part.CreateNewNode(1, 0.0, 0.0, 0.0);

    KRATOS_DEBUG_EXCEPT_EXCEPTION_IS_THROWN(p_node->AddDof(VELOCITY_Y),
        "Error: The Dof-Variable VELOCITY_Y is not in the list of variables");
}

KRATOS_TEST_CASE_IN_SUITE(DofConstructorWtihoutReactionInVariablesList, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& model_part = current_model.CreateModelPart("TestModelPart");

    model_part.AddNodalSolutionStepVariable(VELOCITY);

    model_part.SetBufferSize(1);

    auto p_node = model_part.CreateNewNode(1, 0.0, 0.0, 0.0);

    KRATOS_DEBUG_EXCEPT_EXCEPTION_IS_THROWN(p_node->AddDof(VELOCITY_Y, REACTION_Y),
        "Error: The Reaction-Variable REACTION_Y is not in the list of variables");
}

KRATOS_TEST_CASE_IN_SUITE(DofFixing, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& model_part = current_model.CreateModelPart("TestModelPart");

    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(REACTION);

    model_part.SetBufferSize(1);

    auto p_node = model_part.CreateNewNode(1, 0.0, 0.0, 0.0);

    auto p_dof = p_node->pAddDof(VELOCITY_Y, REACTION_Y);

    // Checking default fixities
    KRATOS_EXPECT_TRUE(p_dof->IsFree());
    KRATOS_EXPECT_FALSE(p_dof->IsFixed());

    // Checking after fixing the dof
    p_dof->FixDof();
    KRATOS_EXPECT_FALSE(p_dof->IsFree());
    KRATOS_EXPECT_TRUE(p_dof->IsFixed());

    // Checking after freeing the dof
    p_dof->FreeDof();
    KRATOS_EXPECT_TRUE(p_dof->IsFree());
    KRATOS_EXPECT_FALSE(p_dof->IsFixed());
}

KRATOS_TEST_CASE_IN_SUITE(DofVariables, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& model_part = current_model.CreateModelPart("TestModelPart");

    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(REACTION);

    model_part.SetBufferSize(1);

    auto p_node = model_part.CreateNewNode(1, 0.0, 0.0, 0.0);

    auto p_dof = p_node->pAddDof(VELOCITY_Y, REACTION_Y);
    auto p_dof_2 = p_node->pAddDof(VELOCITY_Z);

    KRATOS_EXPECT_EQ(VELOCITY_Y, p_dof->GetVariable());
    KRATOS_EXPECT_EQ(REACTION_Y, p_dof->GetReaction());

    KRATOS_EXPECT_FALSE(p_dof_2->HasReaction());
}

KRATOS_TEST_CASE_IN_SUITE(DofEquationId, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& model_part = current_model.CreateModelPart("TestModelPart");

    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(REACTION);

    model_part.SetBufferSize(1);

    auto p_node = model_part.CreateNewNode(1, 0.0, 0.0, 0.0);

    auto p_dof = p_node->pAddDof(VELOCITY_Y, REACTION_Y);

    const std::size_t eq_id= 569;

    p_dof->SetEquationId(eq_id);

    KRATOS_EXPECT_EQ(eq_id, p_dof->EquationId());
}

KRATOS_TEST_CASE_IN_SUITE(DofId, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& model_part = current_model.CreateModelPart("TestModelPart");

    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(REACTION);

    model_part.SetBufferSize(1);

    const std::size_t my_id = 4322;

    auto p_node = model_part.CreateNewNode(my_id, 0.0, 0.0, 0.0);

    auto p_dof = p_node->pAddDof(VELOCITY_Y, REACTION_Y);

    KRATOS_EXPECT_EQ(my_id, p_dof->GetId());

    KRATOS_EXPECT_EQ(p_dof->GetId(), p_dof->Id());
}
}
}  // namespace Kratos.
