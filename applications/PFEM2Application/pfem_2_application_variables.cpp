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

#include "pfem_2_application_variables.h"

namespace Kratos
{

 	KRATOS_CREATE_VARIABLE(double, PRESS_GRADIENT_JUMP)
	KRATOS_CREATE_VARIABLE(double, PRESS_DISCONTINUITY);
	KRATOS_CREATE_VARIABLE(double, INV_LAPLACIAN_ENRICH)
	KRATOS_CREATE_VARIABLE(double, ENRICH_RHS)
	KRATOS_CREATE_VARIABLE(double, G_VALUE)
	KRATOS_CREATE_VARIABLE(double, GRADIENT_DISCONTINUITY)
	KRATOS_CREATE_VARIABLE(double, PREVIOUS_ITERATION_PRESSURE)
	KRATOS_CREATE_VARIABLE(double, FIRST_ITERATION_PRESSURE)
	KRATOS_CREATE_VARIABLE(double, VELOCITY_OVER_ELEM_SIZE)
	KRATOS_CREATE_VARIABLE(double, MEAN_SIZE)
	KRATOS_CREATE_VARIABLE(double, MEAN_VELOCITY_DIFFERENCE)
	KRATOS_CREATE_VARIABLE(double, SPECIFIC_HEAT_CAPACITY_WATER)
	KRATOS_CREATE_VARIABLE(double, SPECIFIC_HEAT_CAPACITY_AIR)
	KRATOS_CREATE_VARIABLE(double, DELTA_TEMPERATURE)
	KRATOS_CREATE_VARIABLE(double, AVAILABLE_AIR_VOLUME)
	KRATOS_CREATE_VARIABLE(double, AVAILABLE_UNBURNED_AIR_VOLUME)
	KRATOS_CREATE_VARIABLE(double, OXYGEN_FRACTION)
	KRATOS_CREATE_VARIABLE(double, CORRECTED_DISTANCE)
	KRATOS_CREATE_VARIABLE(double, SOLID_PRESSURE)
	KRATOS_CREATE_VARIABLE(double, SOLID_YP)
	KRATOS_CREATE_VARIABLE(double, WATER_DISTANCE)
	KRATOS_CREATE_VARIABLE(double, ELASTIC_PRESSURE)
	KRATOS_CREATE_VARIABLE(double, WATER_VOLUME)
	KRATOS_CREATE_VARIABLE(double, VOLUME_CORRECTION)
	KRATOS_CREATE_VARIABLE(double, INLET_VELOCITY)

	KRATOS_CREATE_VARIABLE(bool, USEFUL_ELEMENT_FOR_COMBUSTION)

	KRATOS_CREATE_VARIABLE(Vector, ENRICH_LHS_ROW_3D)
	KRATOS_CREATE_VARIABLE(Vector, WATER_GAUSS_POINT)
    KRATOS_CREATE_VARIABLE(Vector, ELEMENT_MEAN_STRESS)

	//typedef PointerVector< PFEM_Particle, PFEM_Particle*, std::vector<PFEM_Particle*> > ParticlePointerVector;
	KRATOS_CREATE_VARIABLE( ParticlePointerVector , PARTICLE_POINTERS)

	//typedef PointerVector< PFEM_Particle_Fluid, PFEM_Particle_Fluid*, std::vector<PFEM_Particle_Fluid*> > FluidParticlePointerVector;
	KRATOS_CREATE_VARIABLE( FluidParticlePointerVector , FLUID_PARTICLE_POINTERS)

	KRATOS_CREATE_VARIABLE(int, NUMBER_OF_PARTICLES)
	KRATOS_CREATE_VARIABLE(int, NUMBER_OF_PARTICLES_AUX)
	KRATOS_CREATE_VARIABLE(int, NUMBER_OF_WATER_PARTICLES)
	KRATOS_CREATE_VARIABLE(int, NUMBER_OF_FLUID_PARTICLES)
	KRATOS_CREATE_VARIABLE(int, PARTICLE_POINTERS_OFFSET)
	KRATOS_CREATE_VARIABLE(int, WATER_PARTICLE_POINTERS_OFFSET)
	KRATOS_CREATE_VARIABLE(int, USE_PRESS_PROJ)

	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(ENRICH_LHS_ROW)
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(ENRICH_PRESS_PROJ_NEGATIVE)
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(ENRICH_PRESS_PROJ_POSITIVE)
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(SURFACE_NORMAL)
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(SURFACE_COORDINATES)
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(PRESS_PROJ_NO_RO)
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DELTA_VELOCITY)
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(WATER_VELOCITY)
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(WATER_MESH_VELOCITY)
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(PROJECTED_VELOCITY)


	/////////
    KRATOS_CREATE_VARIABLE(double, YP_DISTANCE)
    KRATOS_CREATE_VARIABLE(double, YP_NITSCHE)
    KRATOS_CREATE_VARIABLE(double, PRESSURE_NITSCHE)
    KRATOS_CREATE_VARIABLE(double, NITSCHE_ALPHA)
    KRATOS_CREATE_VARIABLE(double, NITSCHE_DELTA)
    KRATOS_CREATE_VARIABLE(double, TIMEDIS_THETA)
    KRATOS_CREATE_VARIABLE(int, CASE_NUMBER)
    KRATOS_CREATE_VARIABLE(int, ELEMENT_ID)
    KRATOS_CREATE_VARIABLE(double, PROJECTED_DISTANCE)
    KRATOS_CREATE_VARIABLE(double, PROJECTED_DELTA_DISTANCE)
    KRATOS_CREATE_VARIABLE(double, LUMPED_MASS_VALUE)  
    KRATOS_CREATE_VARIABLE(double, LUMPED_MASS_VALUE_NITSCHE) 
    KRATOS_CREATE_VARIABLE(double, R_NODE_DISTANCE)
    KRATOS_CREATE_VARIABLE(double, R_NODE_DISTANCE_NITSCHE)
    KRATOS_CREATE_VARIABLE(double, SCALARPROJECTEDVEL_X)
    KRATOS_CREATE_VARIABLE(double, SCALARPROJECTEDVEL_Y)
    KRATOS_CREATE_VARIABLE(double, SCALARPROJECTEDVEL_Z)
    KRATOS_CREATE_VARIABLE(double, SCALARACCELERATION_X)
    KRATOS_CREATE_VARIABLE(double, SCALARACCELERATION_Y)
    KRATOS_CREATE_VARIABLE(double, SCALARACCELERATION_Z)
    KRATOS_CREATE_VARIABLE(double, PRESSURE_N)
    KRATOS_CREATE_VARIABLE(double, DELTA_DISTANCE_NITSCHE)
    /////////

    /////////
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DELTA_VELOCITY_NITSCHE)
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(VELOCITY_NITSCHE)
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(PROJECTED_VELOCITY_NITSCHE)
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(VELOCITY_BEFORESOLUTION)
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(BODY_FORCE_NITSCHE)  
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(PROJECTED_DELTA_VELOCITY)
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(PROJECTED_DELTA_VELOCITY_NITSCHE)
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(R_NODE_VELOCITY)  
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(R_NODE_VELOCITY_NITSCHE) 
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(MESH_DELTA_VELOCITY)      
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(VELOCITY_PIMPLE)     
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(VELOCITY_N)  
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(BODY_FORCE_N)
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DELTA_ACCELERATION)
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(ACCELERATION_AUX)
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(VELOCITY_AUX)
    /////////


    //////DEBUGVARIABLES/////
    KRATOS_CREATE_VARIABLE(double, ISFOUNDA)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDB)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDONESTEPA)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDONESTEPB)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDRK4A)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDRK4B)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDSUBA)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDSUBB)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDSUB2A)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDSUB2B)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDAFIRST)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDBFIRST)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDASECOND)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDBSECOND)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDATHIRD)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDBTHIRD)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDAFOURTH)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDBFOURTH)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDONESTEPAFIRST)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDONESTEPBFIRST)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDONESTEPASECOND)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDONESTEPBSECOND)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDONESTEPATHIRD)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDONESTEPBTHIRD)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDONESTEPAFOURTH)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDONESTEPBFOURTH)  
    KRATOS_CREATE_VARIABLE(double, ISFOUNDRK4AFIRST)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDRK4BFIRST)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDRK4ASECOND)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDRK4BSECOND)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDRK4ATHIRD)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDRK4BTHIRD)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDRK4AFOURTH)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDRK4BFOURTH)    
    KRATOS_CREATE_VARIABLE(double, ISFOUNDSUBAFIRST)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDSUBBFIRST)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDSUBASECOND)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDSUBBSECOND)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDSUBATHIRD)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDSUBBTHIRD)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDSUBAFOURTH)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDSUBBFOURTH)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDSUB2AFIRST)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDSUB2BFIRST)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDSUB2ASECOND)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDSUB2BSECOND)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDSUB2ATHIRD)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDSUB2BTHIRD)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDSUB2AFOURTH)
    KRATOS_CREATE_VARIABLE(double, ISFOUNDSUB2BFOURTH)   
    

    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(PROJECTED_VELOCITY_FIRST)
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(PROJECTED_VELOCITY_SECOND)
    //////DEBUGVARIABLES/////

} // namespace Kratos
