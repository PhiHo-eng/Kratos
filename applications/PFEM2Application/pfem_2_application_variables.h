//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//


#ifndef KRATOS_PFEM2_APPLICATION_VARIABLES_H_INCLUDED
#define KRATOS_PFEM2_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"

#include "custom_utilities/pfem_particle.h"
#include "custom_utilities/pfem_particle_fluidonly.h"


namespace Kratos
{

    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, PRESS_GRADIENT_JUMP)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, PRESS_DISCONTINUITY)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, INV_LAPLACIAN_ENRICH)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ENRICH_RHS)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, G_VALUE)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, GRADIENT_DISCONTINUITY)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, PREVIOUS_ITERATION_PRESSURE)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, FIRST_ITERATION_PRESSURE)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, VELOCITY_OVER_ELEM_SIZE)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, MEAN_SIZE)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, MEAN_VELOCITY_DIFFERENCE)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, SPECIFIC_HEAT_CAPACITY_WATER)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, SPECIFIC_HEAT_CAPACITY_AIR)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, DELTA_TEMPERATURE)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, AVAILABLE_AIR_VOLUME)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, AVAILABLE_UNBURNED_AIR_VOLUME)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, OXYGEN_FRACTION)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, CORRECTED_DISTANCE)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, SOLID_PRESSURE)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, SOLID_YP)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, WATER_DISTANCE)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ELASTIC_PRESSURE)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, WATER_VOLUME)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, VOLUME_CORRECTION)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, INLET_VELOCITY)

    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,bool, USEFUL_ELEMENT_FOR_COMBUSTION)

    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,Vector, ENRICH_LHS_ROW_3D)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,Vector, WATER_GAUSS_POINT)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,Vector, ELEMENT_MEAN_STRESS)

    typedef PointerVector< PFEM_Particle, PFEM_Particle*, std::vector<PFEM_Particle*> > ParticlePointerVector;
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,ParticlePointerVector , PARTICLE_POINTERS)

    typedef PointerVector< PFEM_Particle_Fluid, PFEM_Particle_Fluid*, std::vector<PFEM_Particle_Fluid*> > FluidParticlePointerVector;
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,FluidParticlePointerVector , FLUID_PARTICLE_POINTERS)

    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,int, NUMBER_OF_PARTICLES)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,int, NUMBER_OF_PARTICLES_AUX)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,int, NUMBER_OF_WATER_PARTICLES)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,int, NUMBER_OF_FLUID_PARTICLES)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,int, PARTICLE_POINTERS_OFFSET)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,int, WATER_PARTICLE_POINTERS_OFFSET)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,int, USE_PRESS_PROJ)

    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(KratosPFEM2Application,ENRICH_LHS_ROW)
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(KratosPFEM2Application,ENRICH_PRESS_PROJ_NEGATIVE)
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(KratosPFEM2Application,ENRICH_PRESS_PROJ_POSITIVE)
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(KratosPFEM2Application,SURFACE_NORMAL)
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(KratosPFEM2Application,SURFACE_COORDINATES)
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(KratosPFEM2Application,PRESS_PROJ_NO_RO)
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(KratosPFEM2Application,DELTA_VELOCITY)
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(KratosPFEM2Application,WATER_VELOCITY)
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(KratosPFEM2Application,WATER_MESH_VELOCITY)
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(KratosPFEM2Application,PROJECTED_VELOCITY)



    /////////
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, YP_DISTANCE)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, YP_NITSCHE)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, PRESSURE_NITSCHE)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, NITSCHE_ALPHA)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, NITSCHE_DELTA)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, TIMEDIS_THETA)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,int, CASE_NUMBER)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,int, ELEMENT_ID)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, PROJECTED_DISTANCE)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, PROJECTED_DELTA_DISTANCE)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, LUMPED_MASS_VALUE) 
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, LUMPED_MASS_VALUE_NITSCHE) 
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, R_NODE_DISTANCE)  
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, R_NODE_DISTANCE_NITSCHE)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, SCALARPROJECTEDVEL_X)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, SCALARPROJECTEDVEL_Y)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, SCALARPROJECTEDVEL_Z)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, SCALARACCELERATION_X)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, SCALARACCELERATION_Y)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, SCALARACCELERATION_Z)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, PRESSURE_N) 
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, DELTA_DISTANCE_NITSCHE) 
    /////////
    

    ////////
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(KratosPFEM2Application,DELTA_VELOCITY_NITSCHE)
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(KratosPFEM2Application,VELOCITY_NITSCHE)
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(KratosPFEM2Application,PROJECTED_VELOCITY_NITSCHE)
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(KratosPFEM2Application,VELOCITY_BEFORESOLUTION)
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(KratosPFEM2Application,BODY_FORCE_NITSCHE)
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(KratosPFEM2Application,PROJECTED_DELTA_VELOCITY)
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(KratosPFEM2Application,PROJECTED_DELTA_VELOCITY_NITSCHE)
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(KratosPFEM2Application,R_NODE_VELOCITY)
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(KratosPFEM2Application,R_NODE_VELOCITY_NITSCHE)
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(KratosPFEM2Application,MESH_DELTA_VELOCITY)
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(KratosPFEM2Application,VELOCITY_PIMPLE)
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(KratosPFEM2Application,VELOCITY_N)
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(KratosPFEM2Application,BODY_FORCE_N)
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(KratosPFEM2Application,DELTA_ACCELERATION)
    /////////

    //////DEBUGVARIABLES/////
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDA)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDB)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDONESTEPA)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDONESTEPB)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDRK4A)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDRK4B)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDSUBA)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDSUBB)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDSUB2A)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDSUB2B)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDAFIRST)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDBFIRST)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDASECOND)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDBSECOND)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDATHIRD)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDBTHIRD)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDAFOURTH)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDBFOURTH)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDONESTEPAFIRST)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDONESTEPBFIRST)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDONESTEPASECOND)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDONESTEPBSECOND)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDONESTEPATHIRD)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDONESTEPBTHIRD)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDONESTEPAFOURTH)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDONESTEPBFOURTH)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDRK4AFIRST)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDRK4BFIRST)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDRK4ASECOND)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDRK4BSECOND)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDRK4ATHIRD)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDRK4BTHIRD)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDRK4AFOURTH)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDRK4BFOURTH)    
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDSUBAFIRST)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDSUBBFIRST)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDSUBASECOND)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDSUBBSECOND)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDSUBATHIRD)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDSUBBTHIRD)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDSUBAFOURTH)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDSUBBFOURTH)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDSUB2AFIRST)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDSUB2BFIRST)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDSUB2ASECOND)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDSUB2BSECOND)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDSUB2ATHIRD)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDSUB2BTHIRD)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDSUB2AFOURTH)
    KRATOS_DEFINE_APPLICATION_VARIABLE(KratosPFEM2Application,double, ISFOUNDSUB2BFOURTH)

    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(KratosPFEM2Application,PROJECTED_VELOCITY_FIRST)
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(KratosPFEM2Application,PROJECTED_VELOCITY_SECOND)
    //////DEBUGVARIABLES/////

}

#endif	// KRATOS_PFEM2_APPLICATION_VARIABLES_H_INCLUDED
