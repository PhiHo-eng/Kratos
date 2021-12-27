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
#include "includes/define.h"
#include "includes/variables.h"
#include "shape_optimization_application.h"

// ==============================================================================

namespace Kratos
{
    // Geometry variables
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(NORMALIZED_SURFACE_NORMAL);

    // Optimization variables
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DF1DX);

    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DF1DX_MAPPED);

    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DC1DX);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DC2DX);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DC3DX);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DC4DX);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DC5DX);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DC6DX);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DC7DX);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DC8DX);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DC9DX);

    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DC1DX_MAPPED);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DC2DX_MAPPED);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DC3DX_MAPPED);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DC4DX_MAPPED);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DC5DX_MAPPED);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DC6DX_MAPPED);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DC7DX_MAPPED);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DC8DX_MAPPED);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DC9DX_MAPPED);

    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(SEARCH_DIRECTION);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(CORRECTION);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(PROJECTION);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(CONTROL_POINT_UPDATE);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(CONTROL_POINT_CHANGE);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(SHAPE_UPDATE);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(SHAPE_CHANGE);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(MESH_CHANGE);

    // For edge damping
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DAMPING_FACTOR);

    // For Mapping
    KRATOS_CREATE_VARIABLE(int,MAPPING_ID);

    // For bead optimization
    KRATOS_CREATE_VARIABLE(double,ALPHA);
    KRATOS_CREATE_VARIABLE(double,ALPHA_MAPPED);
    KRATOS_CREATE_VARIABLE(double,DF1DALPHA);
    KRATOS_CREATE_VARIABLE(double,DF1DALPHA_MAPPED);
    KRATOS_CREATE_VARIABLE(double,DPDALPHA);
    KRATOS_CREATE_VARIABLE(double,DPDALPHA_MAPPED);
    KRATOS_CREATE_VARIABLE(double,DLDALPHA);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(BEAD_DIRECTION);

    // For auxiliary operations
    KRATOS_CREATE_VARIABLE(double,SCALAR_VARIABLE);
    KRATOS_CREATE_VARIABLE(double,SCALAR_VARIABLE_MAPPED);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(VECTOR_VARIABLE);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(VECTOR_VARIABLE_MAPPED);

    // For in plane mapping operations
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(BACKGROUND_COORDINATE);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(BACKGROUND_NORMAL);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(OUT_OF_PLANE_DELTA);

    KRATOS_CREATE_VARIABLE(double, VERTEX_MORPHING_RADIUS);

    // Eof variables

    KratosShapeOptimizationApplication::KratosShapeOptimizationApplication() :
        KratosApplication("ShapeOptimizationApplication")
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

        KRATOS_REGISTER_VARIABLE(VERTEX_MORPHING_RADIUS);
 	}

}  // namespace Kratos.


