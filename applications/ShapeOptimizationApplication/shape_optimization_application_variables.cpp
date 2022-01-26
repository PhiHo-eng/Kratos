//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Reza Najian Asl
//

#include "shape_optimization_application_variables.h"

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

    // For implicit vertex-morphing with Helmholtz PDE
    KRATOS_CREATE_VARIABLE(Matrix, HELMHOLTZ_MASS_MATRIX);
    KRATOS_CREATE_VARIABLE( int,HELMHOLTZ_DIRECTION);
    KRATOS_CREATE_VARIABLE( bool,COMPUTE_CONTROL_POINTS);
    KRATOS_CREATE_VARIABLE( double,HELMHOLTZ_POISSON_RATIO);
    KRATOS_CREATE_VARIABLE( double,HELMHOLTZ_RADIUS);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(HELMHOLTZ_VARS);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(HELMHOLTZ_SOURCE);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(CONTROL_POINT);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(SHAPE);

   // For min-max angles constraints in additive manufacturing 
    KRATOS_CREATE_VARIABLE(double,NODAL_MIN_ANGLE);
    KRATOS_CREATE_VARIABLE(double,NODAL_MAX_ANGLE); 

   // For thickness response
    KRATOS_CREATE_VARIABLE(double,NODAL_THICKNESS); 

}
