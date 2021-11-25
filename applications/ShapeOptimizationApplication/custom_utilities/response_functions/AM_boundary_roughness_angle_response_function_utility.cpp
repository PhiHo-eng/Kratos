// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl
//
// ==============================================================================

// System includes

// External includes

// Project includes
#include "custom_utilities/response_functions/AM_boundary_roughness_angle_response_function_utility.h"
#include "shape_optimization_application.h"
#include "utilities/variable_utils.h"

// ==============================================================================

namespace Kratos
{

AMBoundaryRoughnessAngleResponseFunctionUtility::AMBoundaryRoughnessAngleResponseFunctionUtility(ModelPart& rModelPart, Parameters ResponseSettings)
	: mrModelPart(rModelPart)
{
	const std::size_t domain_size = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];
	KRATOS_ERROR_IF(domain_size != 3) << "AMBoundaryRoughnessAngleResponseFunctionUtility can only be used on 3D geometries!" << std::endl;

	mMainDirection = ResponseSettings["main_direction"].GetVector();
	const double direction_norm = norm_2(mMainDirection);
	KRATOS_ERROR_IF(direction_norm < std::numeric_limits<double>::epsilon()) << "AMBoundaryRoughnessAngleResponseFunctionUtility: 'main_direction' vector norm is 0!" << std::endl;
	mMainDirection /= direction_norm;

	const double max_angle = ResponseSettings["max_angle"].GetDouble();
	mCosMaxAngle = std::cos(max_angle * Globals::Pi / 180);

	const std::string gradient_mode = ResponseSettings["gradient_mode"].GetString();
	if (gradient_mode == "finite_differencing")
	{
		mDelta = ResponseSettings["step_size"].GetDouble();
	}
	else
		KRATOS_ERROR << "Specified gradient_mode '" << gradient_mode << "' not recognized. The only option is: finite_differencing" << std::endl;
}

double AMBoundaryRoughnessAngleResponseFunctionUtility::CalculateValue()
{
	KRATOS_TRY;

	const double value = block_for_each<SumReduction<double>>(mrModelPart.Conditions(), [&](Condition& rCond) {
		return CalculateConditionValue(rCond);});

	mValue = value;

	return mValue;

	KRATOS_CATCH("");
}

void AMBoundaryRoughnessAngleResponseFunctionUtility::CalculateGradient()
{
	KRATOS_TRY;
	// First gradients are initialized
	VariableUtils().SetHistoricalVariableToZero(SHAPE_SENSITIVITY, mrModelPart.Nodes());

	for (auto& cond_i : mrModelPart.Conditions()){

		const double g_i = CalculateConditionValue(cond_i);

		// Compute sensitivities using finite differencing in the three spatial direction
		array_3d gradient(3, 0.0);

		for (auto& node_i : cond_i.GetGeometry()){

			// Apply pertubation in X-direction
			double g_i_after_fd = 0.0;
			node_i.X() += mDelta;
			node_i.X0() += mDelta;
			g_i_after_fd = CalculateConditionValue(cond_i);
			gradient[0] = (g_i_after_fd - g_i) / mDelta;
			node_i.X() -= mDelta;
			node_i.X0() -= mDelta;

			// Apply pertubation in Y-direction
			g_i_after_fd = 0.0;
			node_i.Y() += mDelta;
			node_i.Y0() += mDelta;
			g_i_after_fd = CalculateConditionValue(cond_i);
			gradient[1] = (g_i_after_fd - g_i) / mDelta;
			node_i.Y() -= mDelta;
			node_i.Y0() -= mDelta;

			// Apply pertubation in Z-direction
			g_i_after_fd = 0.0;
			node_i.Z() += mDelta;
			node_i.Z0() += mDelta;
			g_i_after_fd = CalculateConditionValue(cond_i);
			gradient[2] = (g_i_after_fd - g_i) / mDelta;
			node_i.Z() -= mDelta;
			node_i.Z0() -= mDelta;

			// Add to aggregated sensitivities
			noalias(node_i.FastGetSolutionStepValue(SHAPE_SENSITIVITY)) += gradient;
			node_i.FastGetSolutionStepValue(NODAL_MAX_ANGLE) = g_i;
		}
	}

	KRATOS_CATCH("");
}

double AMBoundaryRoughnessAngleResponseFunctionUtility::CalculateConditionValue(const Condition& rFace)
{
	// face normal
	const array_3d local_coords = ZeroVector(3);
	const array_3d face_normal = rFace.GetGeometry().UnitNormal(local_coords);

	double delta_angle = inner_prod(mMainDirection, face_normal) + mCosMaxAngle;
	
	return 1.0/(1.0+std::exp(-20.0*delta_angle));
}

} // namespace Kratos.
