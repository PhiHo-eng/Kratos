//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:
//

// System includes
#include <string>
#include <functional>

// External includes

// Project includes
#include "geometries/geometry.h"
#include "includes/cfd_variables.h"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/kratos_parameters.h"
#include "includes/ublas_interface.h"
#include "response_functions/adjoint_response_function.h"
#include "utilities/element_size_calculator.h"
#include "utilities/geometrical_sensitivity_utility.h"
#include "utilities/parallel_utilities.h"
#include "utilities/variable_utils.h"

// Application includes
#include "custom_utilities/fluid_calculation_utilities.h"

// Include base h
#include "residual_response_function.h"

namespace Kratos
{
template <unsigned int TDim>
ResidualResponseFunction<TDim>::ResidualResponseFunction(
    Parameters Settings,
    ModelPart& rModelPart)
    : mrModelPart(rModelPart)
{
    KRATOS_TRY;

    Parameters default_settings(R"({})");

    Settings.ValidateAndAssignDefaults(default_settings);

    KRATOS_CATCH("");
}

template <unsigned int TDim>
ResidualResponseFunction<TDim>::ResidualResponseFunction(
    ModelPart& rModelPart)
    : mrModelPart(rModelPart)
{
}

template <unsigned int TDim>
void ResidualResponseFunction<TDim>::CalculateGradient(
    const Element& rAdjointElement,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    rResponseGradient.clear();
}

template <unsigned int TDim>
void ResidualResponseFunction<TDim>::CalculateGradient(
    const Condition& rAdjointCondition,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    rResponseGradient.clear();
}

template <unsigned int TDim>
void ResidualResponseFunction<TDim>::CalculateFirstDerivativesGradient(
    const Element& rAdjointElement,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY

    if (rResponseGradient.size() != rResidualGradient.size1()) {
        rResponseGradient.resize(rResidualGradient.size1(), false);
    }

    rResponseGradient.clear();

    const auto& r_geometry = rAdjointElement.GetGeometry();
    const IndexType number_of_nodes = r_geometry.PointsNumber();
    const IndexType block_size = rResponseGradient.size() / number_of_nodes;

    Vector Ws;
    Matrix Ns;
    ShapeFunctionDerivativesArrayType dNdXs;
    this->CalculateGeometryData(rAdjointElement, Ws, Ns, dNdXs);

    const auto& r_properties = rAdjointElement.GetProperties();
    const double density = r_properties[DENSITY];
    const double nu = r_properties[DYNAMIC_VISCOSITY] / density;

    const IndexType number_of_gauss_points = Ws.size();

    MatrixDD velocity_gradient;
    ArrayD velocity, body_force, acceleration, pressure_gradient, gauss_momentum_residual, gauss_momentum_residual_derivative;
    double gauss_continuity_residual;

    double h;
    if (TDim == 2) {
        if (number_of_nodes == 3) {
            h = ElementSizeCalculator<2, 3>::AverageElementSize(r_geometry);
        } else if (number_of_nodes == 4) {
            h = ElementSizeCalculator<2, 4>::AverageElementSize(r_geometry);
        } else {
            KRATOS_ERROR << "Unsupported geometry type having "
                            << number_of_nodes << " nodes in " << TDim << "D.";
        }
    } else if (TDim == 3) {
        if (number_of_nodes == 4) {
            h = ElementSizeCalculator<3, 4>::AverageElementSize(r_geometry);
        } else if (number_of_nodes == 8) {
            h = ElementSizeCalculator<3, 8>::AverageElementSize(r_geometry);
        } else {
            KRATOS_ERROR << "Unsupported geometry type having "
                            << number_of_nodes << " nodes in " << TDim << "D.";
        }
    }

    // calculate the strong form residual
    for (IndexType g = 0; g < number_of_gauss_points; ++g) {
        const double W = Ws[g];
        const Vector& N = row(Ns, g);
        const Matrix& dNdX = dNdXs[g];

        FluidCalculationUtilities::EvaluateInPoint(
            r_geometry, N,
            std::tie(velocity, VELOCITY),
            std::tie(body_force, BODY_FORCE),
            std::tie(acceleration, ACCELERATION));

        FluidCalculationUtilities::EvaluateGradientInPoint(
            r_geometry, dNdX,
            std::tie(velocity_gradient, VELOCITY),
            std::tie(pressure_gradient, PRESSURE));

        const ArrayD& velocity_dot_velocity_gradient = prod(velocity_gradient, velocity);
        const double velocity_magnitude = norm_2(velocity);
        const double inv_velocity_magnitude = (velocity_magnitude > 1e-12) ? 1.0 / velocity_magnitude : 0.0;
        const Vector& convection_operator = prod(dNdX, velocity);

        noalias(gauss_momentum_residual) =
            acceleration + velocity_dot_velocity_gradient +
            pressure_gradient / density - body_force;

        gauss_continuity_residual = 0.0;
        for (IndexType i = 0; i < TDim; ++i) {
            gauss_continuity_residual += velocity_gradient(i, i);
        }

        double tau_u, tau_p;
        CalculateStabilizationParameters(tau_u, tau_p, velocity_magnitude, h, nu, -1.0, rProcessInfo);

        const double gauss_momentum_residual_inner_prod = inner_prod(gauss_momentum_residual, gauss_momentum_residual);

        const double value_1 = std::pow(tau_u, 2) * gauss_momentum_residual_inner_prod;
        const double value_2 = tau_p * gauss_continuity_residual;

        for (IndexType c = 0; c < number_of_nodes; ++c) {
            for (IndexType k = 0; k < TDim; ++k) {
                // compute momentum equation residual derivative
                ArrayD identity_vector = ZeroVector(TDim);
                identity_vector[k] = 1.0;
                noalias(gauss_momentum_residual_derivative) = N[c] * column(velocity_gradient, k) + convection_operator[c] * identity_vector;
                const double gauss_momentum_residual_inner_prod_derivative = 2.0 * inner_prod(gauss_momentum_residual, gauss_momentum_residual_derivative);

                // compute continuity equation residual derivative
                const double gauss_continuity_residual_derivative = dNdX(c, k);

                const double velocity_magnitude_derivative = velocity[k] * N[c] * inv_velocity_magnitude;

                double tau_u_derivative, tau_p_derivative;
                CalculateStabilizationParameterDerivatives(
                    tau_u_derivative, tau_p_derivative, tau_u, tau_p, velocity_magnitude,
                    velocity_magnitude_derivative, h, 0.0, nu, rProcessInfo);

                // derivatives
                double value_1_derivative = 0.0;
                value_1_derivative += gauss_momentum_residual_inner_prod_derivative * std::pow(tau_u, 2);
                value_1_derivative += gauss_momentum_residual_inner_prod * 2.0 * tau_u * tau_u_derivative;

                double value_2_derivative = 0.0;
                value_2_derivative += gauss_continuity_residual_derivative * tau_p;
                value_2_derivative += gauss_continuity_residual * tau_p_derivative;

                double value = 0.0;
                value += 2.0 * value_1 * value_1_derivative;
                value += 2.0 * value_2 * value_2_derivative;

                // add velocity derivatives
                rResponseGradient[c * block_size + k] += value * W;
            }

            noalias(gauss_momentum_residual_derivative) = row(dNdX, c) / density;
            const double gauss_momentum_residual_inner_prod_derivative = 2.0 * inner_prod(gauss_momentum_residual, gauss_momentum_residual_derivative);
            double value_1_derivative = 0.0;
            value_1_derivative += gauss_momentum_residual_inner_prod_derivative * std::pow(tau_u, 2);

            rResponseGradient[c * block_size + TDim] += 2.0 * value_1 * value_1_derivative * W;
        }
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim>
void ResidualResponseFunction<TDim>::CalculateFirstDerivativesGradient(
    const Condition& rAdjointCondition,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    rResponseGradient.clear();
}

template <unsigned int TDim>
void ResidualResponseFunction<TDim>::CalculateSecondDerivativesGradient(
    const Element& rAdjointElement,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY

    if (rResponseGradient.size() != rResidualGradient.size1()) {
        rResponseGradient.resize(rResidualGradient.size1(), false);
    }

    rResponseGradient.clear();

    const auto& r_geometry = rAdjointElement.GetGeometry();
    const IndexType number_of_nodes = r_geometry.PointsNumber();
    const IndexType block_size = rResponseGradient.size() / number_of_nodes;

    Vector Ws;
    Matrix Ns;
    ShapeFunctionDerivativesArrayType dNdXs;
    this->CalculateGeometryData(rAdjointElement, Ws, Ns, dNdXs);

    const auto& r_properties = rAdjointElement.GetProperties();
    const double density = r_properties[DENSITY];
    const double nu = r_properties[DYNAMIC_VISCOSITY] / density;

    double h;
    if (TDim == 2) {
        if (number_of_nodes == 3) {
            h = ElementSizeCalculator<2, 3>::AverageElementSize(r_geometry);
        } else if (number_of_nodes == 4) {
            h = ElementSizeCalculator<2, 4>::AverageElementSize(r_geometry);
        } else {
            KRATOS_ERROR << "Unsupported geometry type having "
                            << number_of_nodes << " nodes in " << TDim << "D.";
        }
    } else if (TDim == 3) {
        if (number_of_nodes == 4) {
            h = ElementSizeCalculator<3, 4>::AverageElementSize(r_geometry);
        } else if (number_of_nodes == 8) {
            h = ElementSizeCalculator<3, 8>::AverageElementSize(r_geometry);
        } else {
            KRATOS_ERROR << "Unsupported geometry type having "
                            << number_of_nodes << " nodes in " << TDim << "D.";
        }
    }

    const IndexType number_of_gauss_points = Ws.size();

    MatrixDD velocity_gradient;
    ArrayD velocity, body_force, acceleration, pressure_gradient, gauss_momentum_residual, gauss_momentum_residual_derivative;
    double gauss_continuity_residual;

    // calculate the strong form residual
    for (IndexType g = 0; g < number_of_gauss_points; ++g) {
        const double W = Ws[g];
        const Vector& N = row(Ns, g);
        const Matrix& dNdX = dNdXs[g];

        FluidCalculationUtilities::EvaluateInPoint(
            r_geometry, N,
            std::tie(velocity, VELOCITY),
            std::tie(body_force, BODY_FORCE),
            std::tie(acceleration, ACCELERATION));

        FluidCalculationUtilities::EvaluateGradientInPoint(
            r_geometry, dNdX,
            std::tie(velocity_gradient, VELOCITY),
            std::tie(pressure_gradient, PRESSURE));

        const ArrayD& velocity_dot_velocity_gradient = prod(velocity_gradient, velocity);

        noalias(gauss_momentum_residual) =
            acceleration + velocity_dot_velocity_gradient +
            pressure_gradient / density - body_force;

        gauss_continuity_residual = 0.0;
        for (IndexType i = 0; i < TDim; ++i) {
            gauss_continuity_residual += velocity_gradient(i, i);
        }

        double tau_u, tau_p;
        CalculateStabilizationParameters(tau_u, tau_p, norm_2(velocity), h, nu, -1.0, rProcessInfo);

        const double gauss_momentum_residual_inner_prod = inner_prod(gauss_momentum_residual, gauss_momentum_residual);

        const double value_1 = std::pow(tau_u, 2) * gauss_momentum_residual_inner_prod;

        for (IndexType c = 0; c < number_of_nodes; ++c) {
            for (IndexType k = 0; k < TDim; ++k) {
                // compute momentum equation residual derivative
                ArrayD identity_vector = ZeroVector(TDim);
                identity_vector[k] = 1.0;
                noalias(gauss_momentum_residual_derivative) = N[c] * identity_vector;
                const double gauss_momentum_residual_inner_prod_derivative = 2.0 * inner_prod(gauss_momentum_residual, gauss_momentum_residual_derivative);

                const double value_1_derivative = std::pow(tau_u, 2) * gauss_momentum_residual_inner_prod_derivative;

                // derivatives
                double value = 0.0;
                value += 2.0 * value_1 * value_1_derivative;

                // add velocity derivatives
                rResponseGradient[c * block_size + k] += value * W;
            }
        }
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim>
void ResidualResponseFunction<TDim>::CalculateSecondDerivativesGradient(
    const Condition& rAdjointCondition,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    rResponseGradient.clear();
}

template <unsigned int TDim>
void ResidualResponseFunction<TDim>::CalculatePartialSensitivity(
    Element& rAdjointElement,
    const Variable<array_1d<double, 3>>& rVariable,
    const Matrix& rSensitivityMatrix,
    Vector& rSensitivityGradient,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY;

    if (rVariable == SHAPE_SENSITIVITY) {
        if (rSensitivityGradient.size() != rSensitivityMatrix.size1()) {
            rSensitivityGradient.resize(rSensitivityMatrix.size1(), false);
        }

        rSensitivityGradient.clear();

        const auto& r_geometry = rAdjointElement.GetGeometry();
        const IndexType number_of_nodes = r_geometry.PointsNumber();

        Vector Ws;
        Matrix Ns;
        ShapeFunctionDerivativesArrayType dNdXs;
        this->CalculateGeometryData(rAdjointElement, Ws, Ns, dNdXs);
        const GeometryData::IntegrationMethod integration_method = rAdjointElement.GetIntegrationMethod();

        const auto& r_properties = rAdjointElement.GetProperties();
        const double density = r_properties[DENSITY];
        const double nu = r_properties[DYNAMIC_VISCOSITY] / density;

        double h;
        std::function<double(unsigned int, unsigned int, const GeometryType&)> element_size_derivative_method;
        if (TDim == 2) {
            if (number_of_nodes == 3) {
                h = ElementSizeCalculator<2, 3>::AverageElementSize(r_geometry);
                element_size_derivative_method = &ElementSizeCalculator<2, 3>::AverageElementSizeDerivative;
            } else if (number_of_nodes == 4) {
                h = ElementSizeCalculator<2, 4>::AverageElementSize(r_geometry);
                element_size_derivative_method = &ElementSizeCalculator<2, 4>::AverageElementSizeDerivative;
            } else {
                KRATOS_ERROR << "Unsupported geometry type having "
                                << number_of_nodes << " nodes in " << TDim << "D.";
            }
        } else if (TDim == 3) {
            if (number_of_nodes == 4) {
                h = ElementSizeCalculator<3, 4>::AverageElementSize(r_geometry);
                element_size_derivative_method = &ElementSizeCalculator<3, 4>::AverageElementSizeDerivative;
            } else if (number_of_nodes == 8) {
                h = ElementSizeCalculator<3, 8>::AverageElementSize(r_geometry);
                element_size_derivative_method = &ElementSizeCalculator<3, 8>::AverageElementSizeDerivative;
            } else {
                KRATOS_ERROR << "Unsupported geometry type having "
                                << number_of_nodes << " nodes in " << TDim << "D.";
            }
        }

        const IndexType number_of_gauss_points = Ws.size();

        MatrixDD velocity_gradient, velocity_gradient_derivative;
        ArrayD velocity, body_force, acceleration, pressure_gradient, pressure_gradient_derivative, gauss_momentum_residual, gauss_momentum_residual_derivative;
        double gauss_continuity_residual;

        // calculate the strong form residual
        for (IndexType g = 0; g < number_of_gauss_points; ++g) {
            const double W = Ws[g];
            const Vector& N = row(Ns, g);
            const Matrix& dNdX = dNdXs[g];

            FluidCalculationUtilities::EvaluateInPoint(
                r_geometry, N,
                std::tie(velocity, VELOCITY),
                std::tie(body_force, BODY_FORCE),
                std::tie(acceleration, ACCELERATION));

            FluidCalculationUtilities::EvaluateGradientInPoint(
                r_geometry, dNdX,
                std::tie(velocity_gradient, VELOCITY),
                std::tie(pressure_gradient, PRESSURE));

            const ArrayD& velocity_dot_velocity_gradient = prod(velocity_gradient, velocity);
            const double velocity_magnitude = norm_2(velocity);

            noalias(gauss_momentum_residual) =
                acceleration + velocity_dot_velocity_gradient +
                pressure_gradient / density - body_force;

            gauss_continuity_residual = 0.0;
            for (IndexType i = 0; i < TDim; ++i) {
                gauss_continuity_residual += velocity_gradient(i, i);
            }

            const double gauss_momentum_residual_inner_prod = inner_prod(gauss_momentum_residual, gauss_momentum_residual);


            double tau_u, tau_p;
            CalculateStabilizationParameters(tau_u, tau_p, velocity_magnitude, h, nu, -1.0, rProcessInfo);

            const double value_1 = std::pow(tau_u, 2) * gauss_momentum_residual_inner_prod;
            const double value_2 = tau_p * gauss_continuity_residual;

            Geometry<Point>::JacobiansType J;
            r_geometry.Jacobian(J, integration_method);
            const auto& DN_De = r_geometry.ShapeFunctionsLocalGradients(integration_method);

            GeometricalSensitivityUtility::ShapeFunctionsGradientType dNdX_deriv;
            const Matrix& rJ = J[g];
            const Matrix& rDN_De = DN_De[g];
            const double inv_detJ = 1.0 / MathUtils<double>::DetMat(rJ);
            GeometricalSensitivityUtility geom_sensitivity(rJ, rDN_De);

            ShapeParameter deriv;
            for (deriv.NodeIndex = 0; deriv.NodeIndex < number_of_nodes; ++deriv.NodeIndex) {
                for (deriv.Direction = 0; deriv.Direction < TDim; ++deriv.Direction) {
                    double detJ_deriv;
                    geom_sensitivity.CalculateSensitivity(deriv, detJ_deriv, dNdX_deriv);
                    const double weight_deriv = detJ_deriv * inv_detJ * W;

                    const double h_derivative = element_size_derivative_method(deriv.NodeIndex, deriv.Direction, r_geometry);
                    double tau_u_derivative, tau_p_derivative;
                    CalculateStabilizationParameterDerivatives(
                        tau_u_derivative, tau_p_derivative, tau_u, tau_p, velocity_magnitude,
                        0.0, h, h_derivative, nu, rProcessInfo);

                    FluidCalculationUtilities::EvaluateGradientInPoint(
                        r_geometry, dNdX_deriv,
                        std::tie(velocity_gradient_derivative, VELOCITY),
                        std::tie(pressure_gradient_derivative, PRESSURE));

                    // compute momentum equation residual derivative
                    noalias(gauss_momentum_residual_derivative) = prod(velocity_gradient_derivative, velocity) + pressure_gradient_derivative / density;
                    const double gauss_momentum_residual_inner_prod_derivative = 2.0 * inner_prod(gauss_momentum_residual, gauss_momentum_residual_derivative);

                    // compute continuity residual derivative
                    double gauss_continuity_residual_derivative = 0.0;
                    for (IndexType i = 0; i < TDim; ++i) {
                        gauss_continuity_residual_derivative += velocity_gradient_derivative(i, i);
                    }

                    // derivatives
                    double value_1_derivative = 0.0;
                    value_1_derivative += gauss_momentum_residual_inner_prod_derivative * std::pow(tau_u, 2);
                    value_1_derivative += gauss_momentum_residual_inner_prod * 2.0 * tau_u * tau_u_derivative;

                    double value_2_derivative = 0.0;
                    value_2_derivative += gauss_continuity_residual_derivative * tau_p;
                    value_2_derivative += gauss_continuity_residual * tau_p_derivative;

                    double value = 0.0;

                    value += 2.0 * value_1 * value_1_derivative * W;
                    value += std::pow(value_1, 2) * weight_deriv;
                    value += 2.0 * value_2 * value_2_derivative * W;
                    value += std::pow(value_2, 2) * weight_deriv;

                    // add velocity derivatives
                    rSensitivityGradient[deriv.NodeIndex * TDim + deriv.Direction] += value;
                }
            }
        }
    } else {
        KRATOS_ERROR << "Unsupported partial sensitivity variable " << rVariable.Name() << " requested.";
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim>
void ResidualResponseFunction<TDim>::CalculatePartialSensitivity(
    Condition& rAdjointCondition,
    const Variable<array_1d<double, 3>>& rVariable,
    const Matrix& rSensitivityMatrix,
    Vector& rSensitivityGradient,
    const ProcessInfo& rProcessInfo)
{
    if (rSensitivityGradient.size() != rSensitivityMatrix.size1())
        rSensitivityGradient.resize(rSensitivityMatrix.size1(), false);

    rSensitivityGradient.clear();
}

template <unsigned int TDim>
double ResidualResponseFunction<TDim>::CalculateValue(ModelPart& rModelPart)
{
    KRATOS_TRY;

    struct TLS {
        Vector Ws;
        Matrix Ns;
        ShapeFunctionDerivativesArrayType dNdXs;
    };

    const double value = block_for_each<SumReduction<double>>(
        mrModelPart.Elements(), TLS(), [&](Element& rElement, TLS& rTLS) {
            const auto& r_geometry = rElement.GetGeometry();
            const IndexType number_of_nodes = r_geometry.PointsNumber();
            this->CalculateGeometryData(rElement, rTLS.Ws, rTLS.Ns, rTLS.dNdXs);

            const auto& r_properties = rElement.GetProperties();
            const double density = r_properties[DENSITY];
            const double nu = r_properties[DYNAMIC_VISCOSITY] / density;
            const IndexType number_of_gauss_points = rTLS.Ws.size();

            MatrixDD velocity_gradient;
            ArrayD velocity, body_force, acceleration, pressure_gradient, gauss_momentum_residual;
            double gauss_continuity_residual;

            double h;
            if (TDim == 2) {
                if (number_of_nodes == 3) {
                    h = ElementSizeCalculator<2, 3>::AverageElementSize(r_geometry);
                } else if (number_of_nodes == 4) {
                    h = ElementSizeCalculator<2, 4>::AverageElementSize(r_geometry);
                } else {
                    KRATOS_ERROR << "Unsupported geometry type having "
                                 << number_of_nodes << " nodes in " << TDim << "D.";
                }
            } else if (TDim == 3) {
                if (number_of_nodes == 4) {
                    h = ElementSizeCalculator<3, 4>::AverageElementSize(r_geometry);
                } else if (number_of_nodes == 8) {
                    h = ElementSizeCalculator<3, 8>::AverageElementSize(r_geometry);
                } else {
                    KRATOS_ERROR << "Unsupported geometry type having "
                                 << number_of_nodes << " nodes in " << TDim << "D.";
                }
            }

            // calculate the strong form residual
            double value;
            for (IndexType g = 0; g < number_of_gauss_points; ++g) {
                const double W = rTLS.Ws[g];
                const Vector& N = row(rTLS.Ns, g);
                const Matrix& dNdX = rTLS.dNdXs[g];

                FluidCalculationUtilities::EvaluateInPoint(
                    r_geometry, N,
                    std::tie(velocity, VELOCITY),
                    std::tie(body_force, BODY_FORCE),
                    std::tie(acceleration, ACCELERATION));

                FluidCalculationUtilities::EvaluateGradientInPoint(
                    r_geometry, dNdX,
                    std::tie(velocity_gradient, VELOCITY),
                    std::tie(pressure_gradient, PRESSURE));

                const ArrayD& velocity_dot_velocity_gradient = prod(velocity_gradient, velocity);

                noalias(gauss_momentum_residual) =
                    acceleration + velocity_dot_velocity_gradient +
                    pressure_gradient / density - body_force;

                gauss_continuity_residual = 0.0;
                for (IndexType i = 0; i < TDim; ++i) {
                    gauss_continuity_residual += velocity_gradient(i, i);
                }

                double tau_u, tau_p;
                CalculateStabilizationParameters(tau_u, tau_p, norm_2(velocity), h, nu, 1.0, mrModelPart.GetProcessInfo());

                // value_1 += std::pow(std::pow(tau_u, 2) * inner_prod(gauss_momentum_residual, gauss_momentum_residual), 2) * W;
                // value_2 += std::pow(tau_p * gauss_continuity_residual, 2) * W;

                value += std::pow(inner_prod(gauss_momentum_residual, gauss_momentum_residual) * gauss_continuity_residual, 2) * W;
            }

            rElement.SetValue(ELEMENT_ERROR, value);

            return value_1 + value_2;
        });

    return mrModelPart.GetCommunicator().GetDataCommunicator().SumAll(value);

    KRATOS_CATCH("");
}

template <unsigned int TDim>
void ResidualResponseFunction<TDim>::CalculateStabilizationParameters(
    double& TauU,
    double& TauP,
    const double VelocityMagnitude,
    const double ElementSize,
    const double KinematicViscosity,
    const double TauDynamicMultiplier,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const double inv_tau_u =
        TauDynamicMultiplier * rCurrentProcessInfo[DYNAMIC_TAU] / rCurrentProcessInfo[DELTA_TIME] +
        2.0 * VelocityMagnitude / ElementSize +
        4.0 * KinematicViscosity / (ElementSize * ElementSize);
    TauU = 1.0 / inv_tau_u;
    TauP = KinematicViscosity + 0.5 * ElementSize * VelocityMagnitude;
}

template <unsigned int TDim>
void ResidualResponseFunction<TDim>::CalculateStabilizationParameterDerivatives(
    double& TauUDerivative,
    double& TauPDerivative,
    const double TauU,
    const double TauP,
    const double VelocityMagnitude,
    const double VelocityMagnitudeDerivative,
    const double ElementSize,
    const double ElementSizeDerivative,
    const double KinematicViscosity,
    const ProcessInfo& rCurrentProcessInfo) const
{
    double inv_tau_u_derivative = 0.0;
    inv_tau_u_derivative += 2.0 * VelocityMagnitudeDerivative / ElementSize;
    inv_tau_u_derivative -= 2.0 * VelocityMagnitude * ElementSizeDerivative / std::pow(ElementSize, 2);
    inv_tau_u_derivative -= 8.0 * KinematicViscosity * ElementSizeDerivative / std::pow(ElementSize, 3);

    TauUDerivative = -1.0 * std::pow(TauU, 2) * inv_tau_u_derivative;

    TauPDerivative = 0.5 * ElementSizeDerivative * VelocityMagnitude;
    TauPDerivative += 0.5 * ElementSize * VelocityMagnitudeDerivative;
}

template <unsigned int TDim>
void ResidualResponseFunction<TDim>::CalculateGeometryData(
    const Element& rElement,
    Vector& rGaussWeights,
    Matrix& rNContainer,
    ShapeFunctionDerivativesArrayType& rDN_DX) const
{
    const GeometryData::IntegrationMethod integration_method =
        rElement.GetIntegrationMethod();
    const GeometryType& r_geometry = rElement.GetGeometry();
    const IndexType number_of_nodes = r_geometry.PointsNumber();
    const IndexType number_of_gauss_points =
        r_geometry.IntegrationPointsNumber(integration_method);

    Vector DetJ;
    r_geometry.ShapeFunctionsIntegrationPointsGradients(rDN_DX, DetJ, integration_method);

    if (rNContainer.size1() != number_of_gauss_points || rNContainer.size2() != number_of_nodes) {
        rNContainer.resize(number_of_gauss_points, number_of_nodes, false);
    }
    rNContainer = r_geometry.ShapeFunctionsValues(integration_method);

    const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
        r_geometry.IntegrationPoints(integration_method);

    if (rGaussWeights.size() != number_of_gauss_points) {
        rGaussWeights.resize(number_of_gauss_points, false);
    }

    for (IndexType g = 0; g < number_of_gauss_points; g++)
        rGaussWeights[g] = DetJ[g] * IntegrationPoints[g].Weight();
}

// template instantiations
template class ResidualResponseFunction<2>;
template class ResidualResponseFunction<3>;

}
