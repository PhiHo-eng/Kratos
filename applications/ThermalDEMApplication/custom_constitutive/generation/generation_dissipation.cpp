//  Kratos Multi-Physics - ThermalDEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

// System includes

// External includes

// Project includes
#include "generation_dissipation.h"

namespace Kratos {
  //-----------------------------------------------------------------------------------------------------------------------
  GenerationDissipation::GenerationDissipation() {}
  GenerationDissipation::~GenerationDissipation() {}

  //------------------------------------------------------------------------------------------------------------
  double GenerationDissipation::ComputeHeatGeneration(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    KRATOS_TRY

    // Check for contact
    if (!particle->mNeighborInContact)
      return 0.0;

    // Add contribution from different sources of energy dissipation
    double heat_generation = 0.0;

    if (r_process_info[GENERATION_SLIDING_OPTION]) {
      heat_generation += ComputeHeatGenerationSlidingFriction(particle);
    }
    if (r_process_info[GENERATION_DAMPING_OPTION]) {
      heat_generation += ComputeHeatGenerationDampingContact(particle);
    }

    // Conversion and partition coefficients
    const double conversion = r_process_info[HEAT_GENERATION_RATIO];
    const double partition  = ComputePartitionCoeff(particle);

    return partition * conversion * heat_generation;

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double GenerationDissipation::ComputeHeatGenerationSlidingFriction(ThermalSphericParticle* particle) {
    KRATOS_TRY

    typename ThermalSphericParticle:: ContactParams contact_params = particle->GetContactParameters();
    const double force_normal     = contact_params.local_force_total[0];
    const double velocity_tangent = contact_params.local_velocity[1];

    if (abs(force_normal)     < std::numeric_limits<double>::epsilon() ||
        abs(velocity_tangent) < std::numeric_limits<double>::epsilon())
      return 0.0;

    const double friction_coeff = particle->GetContactDynamicFrictionCoefficient();
    return friction_coeff * fabs(force_normal * velocity_tangent);

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double GenerationDissipation::ComputeHeatGenerationDampingContact(ThermalSphericParticle* particle) {
    KRATOS_TRY

    typename ThermalSphericParticle:: ContactParams contact_params = particle->GetContactParameters();
    const double force_damp_normal  = fabs(contact_params.local_force_damping[0]);
    const double force_damp_tangent = fabs(contact_params.local_force_damping[1]);
    const double velocity_normal    = fabs(contact_params.local_velocity[0]);
    const double velocity_tangent   = fabs(contact_params.local_velocity[1]);

    return force_damp_normal * velocity_normal + force_damp_tangent * velocity_tangent;

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double GenerationDissipation::ComputePartitionCoeff(ThermalSphericParticle* particle) {
    KRATOS_TRY

    const double k1 = particle->GetParticleConductivity();
    const double k2 = particle->GetNeighborConductivity();
    return k1 / (k1 + k2);

    KRATOS_CATCH("")
  }

} // namespace Kratos
