//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Mohit Tyagi
//

// System includes
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/prism_3d_6.h"
#include "tests/cpp_tests/geometries/test_geometry.h"
#include "utilities/element_size_calculator.h"

namespace Kratos
{
    namespace Testing
    {
        KRATOS_TEST_CASE_IN_SUITE(Prism3D6ElementSize, KratosCoreGeometriesFastSuite)
        {
            // Case 1
            Geometry<NodeType>::PointsArrayType nodes;
            nodes.push_back(NodeType::Pointer(new NodeType(1, 0.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(2, 1.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(3, 0.0, 1.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(4, 0.0, 0.0, 1.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(5, 1.0, 0.0, 1.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(6, 0.0, 1.0, 1.0)));
            const auto prism1 = *GeometryType::Pointer(new Prism3D6<NodeType>(nodes));
            double minimum_size = ElementSizeCalculator<3, 6>::MinimumElementSize(prism1);
            double average_size = ElementSizeCalculator<3, 6>::AverageElementSize(prism1);
            KRATOS_CHECK_NEAR(minimum_size, std::sqrt(0.5), TOLERANCE);
            KRATOS_CHECK_NEAR(average_size, std::pow(0.5, 1.0 / 3.0), TOLERANCE);


            // Case 2
            nodes.clear();
            nodes.push_back(NodeType::Pointer(new NodeType(1, 0.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(2, 1.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(3, 0.0, 1.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(4, 0.0, 0.0, 0.5)));
            nodes.push_back(NodeType::Pointer(new NodeType(5, 1.0, 0.0, 0.5)));
            nodes.push_back(NodeType::Pointer(new NodeType(6, 0.0, 1.0, 0.5)));
            const auto prism2 = *GeometryType::Pointer(new Prism3D6<NodeType>(nodes));
            minimum_size = ElementSizeCalculator<3, 6>::MinimumElementSize(prism2);
            average_size = ElementSizeCalculator<3, 6>::AverageElementSize(prism2);
            KRATOS_CHECK_NEAR(minimum_size, 0.5, TOLERANCE);
            KRATOS_CHECK_NEAR(average_size, std::pow(0.25, 1.0 / 3.0), TOLERANCE); //
        }

    } // namespace Testing

} // namespace Kratos
