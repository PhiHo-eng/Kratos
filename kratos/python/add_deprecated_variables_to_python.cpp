//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "python/add_deprecated_variables_to_python.h"
#include "includes/deprecated_variables.h"

namespace Kratos
{

namespace Python
{

void  AddDeprecatedVariablesToPython(pybind11::module& m)
{
    //for Level Set application:
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,   IS_DUPLICATED )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,   SPLIT_ELEMENT )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,   SPLIT_NODAL )

    //for PFEM fluids application:
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,  IS_JACK_LINK )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,  IMPOSED_PRESSURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,  IMPOSED_VELOCITY_X )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,  IMPOSED_VELOCITY_Y )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,  IMPOSED_VELOCITY_Z )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,  IMPOSED_ANGULAR_VELOCITY_X )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,  IMPOSED_ANGULAR_VELOCITY_Y )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,  IMPOSED_ANGULAR_VELOCITY_Z )

    //For the DEM Application:
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, IMPOSED_VELOCITY_X_VALUE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, IMPOSED_VELOCITY_Y_VALUE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, IMPOSED_VELOCITY_Z_VALUE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, IMPOSED_ANGULAR_VELOCITY_X_VALUE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, IMPOSED_ANGULAR_VELOCITY_Y_VALUE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, IMPOSED_ANGULAR_VELOCITY_Z_VALUE)

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    IS_INLET )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,   IS_INTERFACE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,   IS_VISITED )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, IS_EROSIONABLE )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    IS_STRUCTURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    IS_POROUS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    IS_WATER )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    IS_FLUID )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    IS_BOUNDARY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    IS_FREE_SURFACE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    IS_AIR_EXIT )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    IS_LAGRANGIAN_INLET )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    IS_WATER_ELEMENT )


    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    IS_BURN )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    IS_DRIPPING )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    IS_PERMANENT )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    IS_WALL )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    Ypr ) //var name does not follow standard
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    Yox )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    Yfuel )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    Hfuel )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    Hpr )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    Hpr1 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    Hox )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    RADIATIVE_INTENSITY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    RADIATIVE_INTENSITY_1 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    RADIATIVE_INTENSITY_2 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    RADIATIVE_INTENSITY_3 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    RADIATIVE_INTENSITY_4 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    RADIATIVE_INTENSITY_5 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    RADIATIVE_INTENSITY_6 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    RADIATIVE_INTENSITY_7 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    RADIATIVE_INTENSITY_8 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    RADIATIVE_INTENSITY_9 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    RADIATIVE_INTENSITY_10 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    RADIATIVE_INTENSITY_11 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    RADIATIVE_INTENSITY_12 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    RADIATIVE_INTENSITY_13 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    RADIATIVE_INTENSITY_14 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    RADIATIVE_INTENSITY_15 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    RADIATIVE_INTENSITY_16 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    RADIATIVE_INTENSITY_17 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    RADIATIVE_INTENSITY_18 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    RADIATIVE_INTENSITY_19 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    RADIATIVE_INTENSITY_20 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    RADIATIVE_INTENSITY_21 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    RADIATIVE_INTENSITY_22 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    RADIATIVE_INTENSITY_23 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    RADIATIVE_INTENSITY_24 )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    rhoD )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    xi )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    a )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    b )


    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    IS_SLIP )

    //for Level Set application:
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,    IS_DIVIDED )

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,  xi_c )

	//for Thermo Mechanical Application
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NODE_PROPERTY_ID)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, REF_ID)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PARTICLE_RADIUS)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, POSETIVE_DISTANCE)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NAGATIVE_DISTANCE)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, IS_ESCAPED)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, IS_SOLIDIFIED)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, IS_GRAVITY_FILLING)

}
} // namespace Python.
} // Namespace Kratos
