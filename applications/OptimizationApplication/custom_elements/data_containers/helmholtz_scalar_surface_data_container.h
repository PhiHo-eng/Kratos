//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl
//                   Suneth Warnakulasuriya
//

#pragma once

// System includes
#include <array>

// External includes

// Project includes
#include "geometries/geometry.h"
#include "includes/node.h"
#include "utilities/math_utils.h"

// Application includes
#include "custom_utilities/entity_calculation_utils.h"
#include "optimization_application_variables.h"

namespace Kratos {

template <unsigned int TDim, unsigned int TNumNodes>
class HelmholtzScalarSurfaceDataContainer {
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    using GeometryType = Geometry<Node>;

    static constexpr IndexType NumberOfNodes = TNumNodes;

    static constexpr std::array<const Variable<double>*, 1> TargetVariablesList = {&HELMHOLTZ_SCALAR};

    static constexpr std::array<const Variable<double>*, 1> SourceVariablesList = {&HELMHOLTZ_SCALAR_SOURCE};

    ///@}
    ///@name Public classes
    ///@{

    struct ConstantDataContainer
    {
        const GeometryType& mrGeometry;
        const GeometryData::IntegrationMethod& mrIntegrationMethod;
        BoundedMatrix<double, 3, 3> mTangentProjectionMatrix;

        ConstantDataContainer(
            const GeometryType& rGeometry,
            const GeometryData::IntegrationMethod& rIntegrationMethod)
            : mrGeometry(rGeometry),
              mrIntegrationMethod(rIntegrationMethod)
        {
        }
    };

    ///@}
    ///@name Life cycle
    ///@{

    HelmholtzScalarSurfaceDataContainer(const Geometry<Node>& rGeometry)
        : mpSolidGeometry(EntityCalculationUtils::CreateSolidGeometry(rGeometry))
    {
    }

    void CalculateConstants(ConstantDataContainer& rConstantData) const
    {
        const auto& integration_points = rConstantData.mrGeometry.IntegrationPoints(rConstantData.mrIntegrationMethod);

        array_1d<double, 3> normal = ZeroVector(3);

        for (IndexType point_number = 0;
             point_number < integration_points.size(); ++point_number) {
            normal += rConstantData.mrGeometry.UnitNormal(point_number, rConstantData.mrIntegrationMethod);
        }

        normal /= integration_points.size();
        normal /= MathUtils<double>::Norm3(normal);

        const BoundedMatrix<double, 3, 3>& id_matrix = IdentityMatrix(3, 3);

        noalias(rConstantData.mTangentProjectionMatrix) = id_matrix - outer_prod(normal, normal);
    }

    void CalculateShapeFunctionDerivatives(
        Matrix& rdNdX,
        const IndexType IntegrationPoint,
        const ConstantDataContainer& rConstantData) const
    {
        Matrix DN_DX;
        EntityCalculationUtils::CalculateSurfaceElementShapeDerivatives(DN_DX, *mpSolidGeometry, rConstantData.mrGeometry, rConstantData.mrIntegrationMethod, IntegrationPoint);
        rdNdX = prod(DN_DX, rConstantData.mTangentProjectionMatrix);
    }

    ///@}

private:
    ///@name Private member variables
    ///@{

    const Geometry<Node>::Pointer mpSolidGeometry;

    ///@}
};

} // namespace Kratos