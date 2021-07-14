// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:         BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if !defined(KRATOS_CONVECTION_DIFUSSION_APPLICATION_VARIABLES_H_INCLUDED)
#define KRATOS_CONVECTION_DIFUSSION_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/global_pointer_variables.h"

namespace Kratos
{

// Variables definition
KRATOS_DEFINE_APPLICATION_VARIABLE( CONVECTION_DIFFUSION_APPLICATION, double, AUX_FLUX)
KRATOS_DEFINE_APPLICATION_VARIABLE( CONVECTION_DIFFUSION_APPLICATION, double, AUX_TEMPERATURE)
KRATOS_DEFINE_APPLICATION_VARIABLE( CONVECTION_DIFFUSION_APPLICATION, double, BFECC_ERROR)
KRATOS_DEFINE_APPLICATION_VARIABLE( CONVECTION_DIFFUSION_APPLICATION, double, BFECC_ERROR_1)
KRATOS_DEFINE_APPLICATION_VARIABLE( CONVECTION_DIFFUSION_APPLICATION, double, DELTA_SCALAR1)
KRATOS_DEFINE_APPLICATION_VARIABLE( CONVECTION_DIFFUSION_APPLICATION, double, MEAN_SIZE)
KRATOS_DEFINE_APPLICATION_VARIABLE( CONVECTION_DIFFUSION_APPLICATION, double, MEAN_VEL_OVER_ELEM_SIZE)
KRATOS_DEFINE_APPLICATION_VARIABLE( CONVECTION_DIFFUSION_APPLICATION, double, MELT_TEMPERATURE_1)
KRATOS_DEFINE_APPLICATION_VARIABLE( CONVECTION_DIFFUSION_APPLICATION, double, MELT_TEMPERATURE_2)
KRATOS_DEFINE_APPLICATION_VARIABLE( CONVECTION_DIFFUSION_APPLICATION, double, PROJECTED_SCALAR1)
KRATOS_DEFINE_APPLICATION_VARIABLE( CONVECTION_DIFFUSION_APPLICATION, double, TRANSFER_COEFFICIENT)
KRATOS_DEFINE_APPLICATION_VARIABLE( CONVECTION_DIFFUSION_APPLICATION, double, ADJOINT_HEAT_TRANSFER)
KRATOS_DEFINE_APPLICATION_VARIABLE( CONVECTION_DIFFUSION_APPLICATION, double, SCALAR_PROJECTION)

KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( CONVECTION_DIFFUSION_APPLICATION, CONVECTION_VELOCITY)

////////
KRATOS_DEFINE_APPLICATION_VARIABLE(CONVECTION_DIFFUSION_APPLICATION,double, SCALARPROJECTEDVEL_X)
KRATOS_DEFINE_APPLICATION_VARIABLE(CONVECTION_DIFFUSION_APPLICATION,double, SCALARPROJECTEDVEL_Y)
KRATOS_DEFINE_APPLICATION_VARIABLE(CONVECTION_DIFFUSION_APPLICATION,double, SCALARPROJECTEDVEL_Z)
KRATOS_DEFINE_APPLICATION_VARIABLE(CONVECTION_DIFFUSION_APPLICATION,double, SCALARACCELERATION_X)
KRATOS_DEFINE_APPLICATION_VARIABLE(CONVECTION_DIFFUSION_APPLICATION,double, SCALARACCELERATION_Y)
KRATOS_DEFINE_APPLICATION_VARIABLE(CONVECTION_DIFFUSION_APPLICATION,double, SCALARACCELERATION_Z)
KRATOS_DEFINE_APPLICATION_VARIABLE(CONVECTION_DIFFUSION_APPLICATION,double, LUMPED_MASS_VALUE)
////////    


////////    
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(CONVECTION_DIFFUSION_APPLICATION,DELTA_ACCELERATION)
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(CONVECTION_DIFFUSION_APPLICATION,ACCELERATION_AUX)
////////    

}

#endif /* KRATOS_CONVECTION_DIFUSSION_APPLICATION_VARIABLES_H_INCLUDED */
