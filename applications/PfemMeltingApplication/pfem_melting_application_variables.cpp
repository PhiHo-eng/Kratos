// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:         BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Julio Marti
//

#include "pfem_melting_application_variables.h"

namespace Kratos
{

KRATOS_CREATE_VARIABLE(double, ACTIVATION_ENERGY)
KRATOS_CREATE_VARIABLE(double, ARRHENIUS_COEFFICIENT)
KRATOS_CREATE_VARIABLE(double, RADIOUS)
KRATOS_CREATE_VARIABLE(double, HEAT_OF_VAPORIZATION)

KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(INITIAL_POSITION)

}
