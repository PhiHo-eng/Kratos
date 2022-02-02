//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Javier San Mauro Saiz
//                   Joaquin Irazabal Gonzalez
//

// Application includes
#include "custom_constitutive/joint_cohesion_driven_2D_law.hpp"

namespace Kratos
{

void JointCohesionDriven2DLaw::GetLawFeatures(Features& rFeatures)
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

void JointCohesionDriven2DLaw::ComputeEquivalentStrain(ConstitutiveLawVariables& rVariables,
                                                       Parameters& rValues)
{
	rVariables.EquivalentStrain = 0.0;

	const Vector& StrainVector = rValues.GetStrainVector();

	if (mStateVariable == 0.0) // Unbroken joint
    {
		// Triangular broken law

        if( rValues.GetOptions().Is(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY) ) // No contact between interfaces
        {
			double tau = rVariables.YoungModulus * StrainVector[0];
			double sigma = rVariables.YoungModulus * StrainVector[1];

			if (sigma > 0.0)
			{
		        double broken_limit = -rVariables.FrictionCoefficient * sigma + rVariables.Cohesion;

    		    if (sigma > rVariables.Cohesion) rVariables.EquivalentStrain = 1.0;
	    	    if (abs (tau) > broken_limit) rVariables.EquivalentStrain = 1.0;
			}
		}
	}
}

//----------------------------------------------------------------------------------------

void JointCohesionDriven2DLaw::ComputeConstitutiveMatrix(Matrix& rConstitutiveMatrix,
                                                         ConstitutiveLawVariables& rVariables,
                                                         Parameters& rValues)
{
    const Vector& StrainVector = rValues.GetStrainVector();

	double broken_YieldStress = 1e-3 * rVariables.YoungModulus;

    if( rValues.GetOptions().Is(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY) ) // No contact between interfaces
    {
        // Tensile constitutive matrix
        if (mStateVariable == 0.0) // Unbroken joint
        {
			rConstitutiveMatrix(0,0) = rVariables.YoungModulus;
			rConstitutiveMatrix(1,1) = rConstitutiveMatrix(0,0);
			rConstitutiveMatrix(0,1) = 0.0;
			rConstitutiveMatrix(1,0) = 0.0;
		}

		else // Broken joint
		{
			rConstitutiveMatrix(0,0) = broken_YieldStress;
			rConstitutiveMatrix(1,1) = rConstitutiveMatrix(0,0);
			rConstitutiveMatrix(0,1) = 0.0;
			rConstitutiveMatrix(1,0) = 0.0;
		}
    }

    else // Contact between interfaces
    {
        // Compresive constitutive matrix
        if (mStateVariable == 0.0) // Unbroken joint
		{
			rConstitutiveMatrix(0,0) = rVariables.YoungModulus;
			rConstitutiveMatrix(1,1) = rConstitutiveMatrix(0,0);
			rConstitutiveMatrix(0,1) = 0.0;
			rConstitutiveMatrix(1,0) = 0.0;
		}

		else // Broken joint
		{
			double shear_modulus_stress = fabs(StrainVector[0] / (2.0 * (1.0 + rVariables.PoissonCoefficient)));
			double friction_modulus_stress = fabs(rVariables.FrictionCoefficient * StrainVector[1]);
			double shear_modulus = rVariables.YoungModulus / (2.0 * (1.0 + rVariables.PoissonCoefficient));

			if (shear_modulus_stress > friction_modulus_stress)
			{
				rConstitutiveMatrix(0,0) = shear_modulus;
				rConstitutiveMatrix(1,1) = rVariables.YoungModulus;
				rConstitutiveMatrix(1,0) = 0.0;

				const double eps = std::numeric_limits<double>::epsilon();
				if(StrainVector[0] > eps)
				{
					rConstitutiveMatrix(0,1) =-rVariables.YoungModulus * rVariables.FrictionCoefficient;
				}
				else if(StrainVector[0] < -eps)
				{
					rConstitutiveMatrix(0,1) = rVariables.YoungModulus * rVariables.FrictionCoefficient;
				}
				else
				{
					rConstitutiveMatrix(0,1) = 0.0;
				}
			}

			if (shear_modulus_stress <= friction_modulus_stress)
			{
				rConstitutiveMatrix(0,0) = shear_modulus;
				rConstitutiveMatrix(1,1) = rVariables.YoungModulus;
				rConstitutiveMatrix(0,1) = 0.0;
				rConstitutiveMatrix(1,0) = 0.0;
			}
		}
    }
}

//----------------------------------------------------------------------------------------

void JointCohesionDriven2DLaw::ComputeStressVector(Vector& rStressVector,
                                                   ConstitutiveLawVariables& rVariables,
                                                   Parameters& rValues)
{
    const Vector& StrainVector = rValues.GetStrainVector();

	double broken_YieldStress = 1e-3 * rVariables.YoungModulus;

    if( rValues.GetOptions().Is(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY) ) // No contact between interfaces
    {
		// Tensile stress
        if (mStateVariable == 0.0) // Unbroken joint
		{
			rStressVector[0] = rVariables.YoungModulus * StrainVector[0];
			rStressVector[1] = rVariables.YoungModulus * StrainVector[1];
		}

		else // Broken joint
		{
			rStressVector[0] = broken_YieldStress * StrainVector[0];
			rStressVector[1] = broken_YieldStress * StrainVector[1];
		}
    }

    else // Contact between interfaces
    {
        // Note: StrainVector[1] < 0.0, rStressVector[1] < 0.0 -> Compresive stress
		if (mStateVariable == 0.0) // Unbroken joint
		{
			rStressVector[0] = rVariables.YoungModulus * StrainVector[0];
			rStressVector[1] = rVariables.YoungModulus * StrainVector[1];
		}

        else // Broken joint
		{
			const double shear_modulus = rVariables.YoungModulus / (2.0 * (1.0 + rVariables.PoissonCoefficient));
			double friction_stress = fabs(shear_modulus * StrainVector[0]);
			double max_friction_stress = fabs(rVariables.FrictionCoefficient * rVariables.YoungModulus * StrainVector[1]);
			if (friction_stress > max_friction_stress) friction_stress = max_friction_stress;

			const double eps = std::numeric_limits<double>::epsilon();
			if(StrainVector[0] > eps)
			{
				rStressVector[0] = + friction_stress;
			}
			else if(StrainVector[0] < -eps)
			{
				rStressVector[0] = - friction_stress;
			}
			else
			{
				rStressVector[0] = 0.0;
			}
			rStressVector[1] = rVariables.YoungModulus * StrainVector[1];
		}
    }
    // Add Uplift Pressure
    rStressVector[1] -= mUpliftPressure;
}

} // Namespace Kratos
