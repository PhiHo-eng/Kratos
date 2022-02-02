//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//

// Application includes
#include "custom_constitutive/joint_bilinear_cohesive_2D_law.hpp"

namespace Kratos
{

void JointBilinearCohesive2DLaw::GetLawFeatures(Features& rFeatures)
{
    //Set the type of law
	rFeatures.mOptions.Set( PLANE_STRAIN_LAW );
	rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
	rFeatures.mOptions.Set( ISOTROPIC );

	//Set strain measure required by the consitutive law
	rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
	//rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

	//Set the spacedimension
	rFeatures.mSpaceDimension = 2;

	//Set the strain size
	rFeatures.mStrainSize = 2;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void JointBilinearCohesive2DLaw::ComputeEquivalentStrain(ConstitutiveLawVariables& rVariables,
                                                         Parameters& rValues)
{
    const Vector& StrainVector = rValues.GetStrainVector();

    if( rValues.GetOptions().Is(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY) ) // No contact between interfaces
    {
		rVariables.EquivalentStrain = std::abs(StrainVector[1]) / rVariables.CriticalDisplacement;
	}
	else // Contact between interfaces
	{
		rVariables.EquivalentStrain = std::abs(StrainVector[0]) / rVariables.CriticalDisplacementTangent;
	}
}

//----------------------------------------------------------------------------------------

void JointBilinearCohesive2DLaw::ComputeConstitutiveMatrix(Matrix& rConstitutiveMatrix,
                                                           ConstitutiveLawVariables& rVariables,
                                                           Parameters& rValues)
{
    const Vector& StrainVector = rValues.GetStrainVector();

    if( rValues.GetOptions().Is(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY) ) // No contact between interfaces
    {
        double aux_param = (1.0 - mStateVariable) / ((1.0 - rVariables.DamageThreshold) * mStateVariable * rVariables.CriticalDisplacement);

        rConstitutiveMatrix(0,0) = rVariables.YieldStress * aux_param;
        rConstitutiveMatrix(1,1) = rConstitutiveMatrix(0,0);
        rConstitutiveMatrix(0,1) = 0.0;
        rConstitutiveMatrix(1,0) = 0.0;
    }
    else // Contact between interfaces
    {
        double aux_param_tangent = (1.0 - mStateVariable) / ((1.0 - rVariables.DamageThreshold) * mStateVariable * rVariables.CriticalDisplacementTangent);
        double aux_param_no_damage = 1.0 / (rVariables.DamageThreshold * rVariables.CriticalDisplacement);

		// Note: StrainVector[1] < 0.0, rStressVector[1] < 0.0 -> Compresive stress
		double actual_friction_stress = fabs(StrainVector[0] * rVariables.YieldStressTangent * aux_param_tangent);
		double max_friction_stress = fabs(rVariables.FrictionCoefficient * rVariables.YoungModulus * aux_param_no_damage * StrainVector[1]);

		if (actual_friction_stress > max_friction_stress)
		{
			rConstitutiveMatrix(0,0) = rVariables.YieldStressTangent * aux_param_tangent;
			rConstitutiveMatrix(1,1) = rVariables.YoungModulus * aux_param_no_damage;
			rConstitutiveMatrix(0,1) = 0.0;
			rConstitutiveMatrix(1,0) = 0.0;
		}
		else
		{
			double friction;
		    if (max_friction_stress < fabs(StrainVector[0] * rVariables.YieldStressTangent * aux_param_no_damage))
			{
                friction = rVariables.YoungModulus * aux_param_no_damage * rVariables.FrictionCoefficient;
			}
            else
			{
				friction = rVariables.YieldStressTangent * aux_param_no_damage * fabs(StrainVector[0] / StrainVector[1]);
			}

			rConstitutiveMatrix(0,0) = rVariables.YieldStressTangent * aux_param_tangent;
			rConstitutiveMatrix(1,1) = rVariables.YoungModulus * aux_param_no_damage;

			const double eps = std::numeric_limits<double>::epsilon();
			if(StrainVector[0] > eps)
			{
				rConstitutiveMatrix(0,1) = -friction;
			}
			else if(StrainVector[0] < -eps)
			{
				rConstitutiveMatrix(0,1) = +friction;
			}
			else
			{
				rConstitutiveMatrix(0,1) = 0.0;
			}
			rConstitutiveMatrix(1,0) = 0.0;
		}
    }
}

//----------------------------------------------------------------------------------------

void JointBilinearCohesive2DLaw::ComputeStressVector(Vector& rStressVector,
                                                     ConstitutiveLawVariables& rVariables,
                                                     Parameters& rValues)
{
    const Vector& StrainVector = rValues.GetStrainVector();

    if( rValues.GetOptions().Is(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY) ) // No contact between interfaces
    {
        double aux_param = (1.0 - mStateVariable) / ((1.0 - rVariables.DamageThreshold) * mStateVariable * rVariables.CriticalDisplacement);

        rStressVector[0] = rVariables.YieldStress * aux_param * StrainVector[0];
        rStressVector[1] = rVariables.YieldStress * aux_param * StrainVector[1];
    }

    else // Contact between interfaces
    {
        double aux_param_tangent = (1.0 - mStateVariable) / ((1.0 - rVariables.DamageThreshold) * mStateVariable * rVariables.CriticalDisplacementTangent);
        double aux_param_no_damage = 1.0 / (rVariables.DamageThreshold * rVariables.CriticalDisplacement);

		rStressVector[1] = rVariables.YoungModulus * aux_param_no_damage * StrainVector[1];

        // Note: StrainVector[1] < 0.0, rStressVector[1] < 0.0 -> Compresive stress
		double actual_friction_stress = fabs(StrainVector[0] * rVariables.YieldStressTangent * aux_param_tangent);
		double max_friction_stress = fabs(rVariables.FrictionCoefficient * rStressVector[1]);

        if (actual_friction_stress > max_friction_stress)
		{
			rStressVector[0] = StrainVector[0] * rVariables.YieldStressTangent * aux_param_tangent;
		}
		else
		{
		    if (max_friction_stress < fabs(StrainVector[0] * rVariables.YieldStressTangent * aux_param_no_damage))
			{
			    double friction_stress;
		        if (max_friction_stress < fabs(StrainVector[0] * rVariables.YieldStressTangent * aux_param_no_damage))
			    {
                    friction_stress = rVariables.FrictionCoefficient * rStressVector[1];
			    }
                else
			    {
				    friction_stress = rVariables.YieldStressTangent * aux_param_no_damage * fabs(StrainVector[0] / StrainVector[1]) * StrainVector[1];
			    }
				const double eps = std::numeric_limits<double>::epsilon();
				if(StrainVector[0] > eps)
				{
					rStressVector[0] = rVariables.YieldStressTangent * aux_param_tangent * StrainVector[0] - friction_stress;
				}
				else if(StrainVector[0] < -eps)
				{
					rStressVector[0] = rVariables.YieldStressTangent * aux_param_tangent * StrainVector[0] + friction_stress;
				}
				else
				{
					rStressVector[0] = 0.0;
				}
			}
			else
			{
				rStressVector[0] = StrainVector[0] * rVariables.YieldStressTangent * aux_param_no_damage;
			}
		}
    }
    // Add Uplift Pressure
    rStressVector[1] -= mUpliftPressure;
}

} // Namespace Kratos
