// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//                   Geiser Armin, https://github.com/armingeiser
//
// ==============================================================================

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_27.h"
#include "shape_optimization_application.h"
#include "shape_optimization_application_variables.h"

// ==============================================================================

namespace Kratos
{
    KratosShapeOptimizationApplication::KratosShapeOptimizationApplication() :
        KratosApplication("ShapeOptimizationApplication"),
            mHelmholtz2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
            mHelmholtz3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
            mHelmholtz3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8<Node<3> >(Element::GeometryType::PointsArrayType(8)))),
            mHelmholtz3D27N(0, Element::GeometryType::Pointer(new Hexahedra3D27<Node<3> >(Element::GeometryType::PointsArrayType(27)))),
            mHelmholtzVec2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
            mHelmholtzVec3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
            mHelmholtzVec3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8<Node<3> >(Element::GeometryType::PointsArrayType(8)))),
            mHelmholtzVec3D27N(0, Element::GeometryType::Pointer(new Hexahedra3D27<Node<3> >(Element::GeometryType::PointsArrayType(27))))
    {}

 	void KratosShapeOptimizationApplication::Register()
 	{
        KRATOS_INFO("") << "    KRATOS   __| |  |   \\   _ \\ __|\n"
                        << "           \\__ \\ __ |  _ \\  __/ _|\n"
                        << "           ____/_| _|_/  _\\_|  ___| OPTIMIZATION\n"
                        << "Initializing KratosShapeOptimizationApplication..." << std::endl;

        // Register variables

        // Geometry variables
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(NORMALIZED_SURFACE_NORMAL);

        // Optimization variables
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DF1DX);

        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DF1DX_MAPPED);

        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DC1DX);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DC2DX);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DC3DX);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DC4DX);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DC5DX);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DC6DX);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DC7DX);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DC8DX);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DC9DX);

        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DC1DX_MAPPED);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DC2DX_MAPPED);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DC3DX_MAPPED);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DC4DX_MAPPED);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DC5DX_MAPPED);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DC6DX_MAPPED);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DC7DX_MAPPED);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DC8DX_MAPPED);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DC9DX_MAPPED);

        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(SEARCH_DIRECTION);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(CORRECTION);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(CONTROL_POINT_UPDATE);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(CONTROL_POINT_CHANGE);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(SHAPE_UPDATE);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(SHAPE_CHANGE);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(MESH_CHANGE);

        // For edge damping
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DAMPING_FACTOR);

        // For mapping
        KRATOS_REGISTER_VARIABLE(MAPPING_ID);

        // For bead optimization
        KRATOS_REGISTER_VARIABLE(ALPHA);
        KRATOS_REGISTER_VARIABLE(ALPHA_MAPPED);
        KRATOS_REGISTER_VARIABLE(DF1DALPHA);
        KRATOS_REGISTER_VARIABLE(DF1DALPHA_MAPPED);
        KRATOS_REGISTER_VARIABLE(DPDALPHA);
        KRATOS_REGISTER_VARIABLE(DPDALPHA_MAPPED);
        KRATOS_REGISTER_VARIABLE(DLDALPHA);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(BEAD_DIRECTION);

        // For auxiliary operations
        KRATOS_REGISTER_VARIABLE(SCALAR_VARIABLE);
        KRATOS_REGISTER_VARIABLE(SCALAR_VARIABLE_MAPPED);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(VECTOR_VARIABLE);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(VECTOR_VARIABLE_MAPPED);

        // For in plane mapping operations
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(BACKGROUND_COORDINATE);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(BACKGROUND_NORMAL);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(OUT_OF_PLANE_DELTA);

        // For implicit vertex-morphing with Helmholtz PDE
        KRATOS_REGISTER_VARIABLE(HELMHOLTZ_DIRECTION)
        KRATOS_REGISTER_VARIABLE(HELMHOLTZ_POISSON_RATIO)
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(HELMHOLTZ_VARS)
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(HELMHOLTZ_SOURCE)        

        // Shape optimization elements

        KRATOS_REGISTER_ELEMENT("HelmholtzElement2D3N", mHelmholtz2D3N);
        KRATOS_REGISTER_ELEMENT("HelmholtzElement3D4N", mHelmholtz3D4N);
        KRATOS_REGISTER_ELEMENT("HelmholtzElement3D8N", mHelmholtz3D8N);
        KRATOS_REGISTER_ELEMENT("HelmholtzElement3D27N", mHelmholtz3D27N);

        KRATOS_REGISTER_ELEMENT("HelmholtzVecElement2D3N", mHelmholtzVec2D3N);
        KRATOS_REGISTER_ELEMENT("HelmholtzVecElement3D4N", mHelmholtzVec3D4N);
        KRATOS_REGISTER_ELEMENT("HelmholtzVecElement3D8N", mHelmholtzVec3D8N);
        KRATOS_REGISTER_ELEMENT("HelmholtzVecElement3D27N",mHelmholtzVec3D27N); 
 	}

}  // namespace Kratos.


