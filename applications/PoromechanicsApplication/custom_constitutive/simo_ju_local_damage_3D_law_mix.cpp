//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

// Application includes
#include "custom_constitutive/simo_ju_local_damage_3D_law_mix.hpp"

namespace Kratos
{

//Default Constructor
SimoJuLocalDamage3DLawMix::SimoJuLocalDamage3DLawMix() : SimoJuLocalDamage3DLaw() {}

//----------------------------------------------------------------------------------------

//Second Constructor
SimoJuLocalDamage3DLawMix::SimoJuLocalDamage3DLawMix(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
    : SimoJuLocalDamage3DLaw(pFlowRule, pYieldCriterion, pHardeningLaw) {}

//----------------------------------------------------------------------------------------

//Copy Constructor
SimoJuLocalDamage3DLawMix::SimoJuLocalDamage3DLawMix(const SimoJuLocalDamage3DLawMix& rOther) : SimoJuLocalDamage3DLaw(rOther) {}

//----------------------------------------------------------------------------------------

//Destructor
SimoJuLocalDamage3DLawMix::~SimoJuLocalDamage3DLawMix() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

ConstitutiveLaw::Pointer SimoJuLocalDamage3DLawMix::Clone() const
{
    SimoJuLocalDamage3DLawMix::Pointer p_clone(new SimoJuLocalDamage3DLawMix(*this));
    return p_clone;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


void LinearElasticPlastic3DLaw::CalculateLinearElasticMatrix( Matrix& rLinearElasticMatrix,const double& YoungModulus,const double& PoissonCoefficient )
{
    rLinearElasticMatrix.clear();

    // 3D linear elastic constitutive matrix
    // LinearElasticMatrix of glass fiber
    Matrix LinearElasticMatrix_glass(6,6);
    noalias(LinearElasticMatrix_glass) = ZeroMatrix(6,6);
    LinearElasticMatrix_glass ( 0 , 0 ) = (YoungModulus*(1.0-PoissonCoefficient)/((1.0+PoissonCoefficient)*(1.0-2.0*PoissonCoefficient)));
    LinearElasticMatrix_glass ( 1 , 1 ) = LinearElasticMatrix_glass ( 0 , 0 );
    LinearElasticMatrix_glass ( 2 , 2 ) = LinearElasticMatrix_glass ( 0 , 0 );

    LinearElasticMatrix_glass ( 3 , 3 ) = LinearElasticMatrix_glass ( 0 , 0 )*(1.0-2.0*PoissonCoefficient)/(2.0*(1.0-PoissonCoefficient));
    LinearElasticMatrix_glass ( 4 , 4 ) = LinearElasticMatrix_glass ( 3 , 3 );
    LinearElasticMatrix_glass ( 5 , 5 ) = LinearElasticMatrix_glass ( 3 , 3 );

    LinearElasticMatrix_glass ( 0 , 1 ) = LinearElasticMatrix_glass ( 0 , 0 )*PoissonCoefficient/(1.0-PoissonCoefficient);
    LinearElasticMatrix_glass ( 1 , 0 ) = LinearElasticMatrix_glass ( 0 , 1 );

    LinearElasticMatrix_glass ( 0 , 2 ) = LinearElasticMatrix_glass ( 0 , 1 );
    LinearElasticMatrix_glass ( 2 , 0 ) = LinearElasticMatrix_glass ( 0 , 1 );

    LinearElasticMatrix_glass ( 1 , 2 ) = LinearElasticMatrix_glass ( 0 , 1 );
    LinearElasticMatrix_glass ( 2 , 1 ) = LinearElasticMatrix_glass ( 0 , 1 );

    // LinearElasticMatrix of Epoxy
    const double young_epoxy = 3.17e9;
    const double poisson_epoxy = 0.38;
    Matrix LinearElasticMatrix_epoxy(6,6);
    noalias(LinearElasticMatrix_epoxy) = ZeroMatrix(6,6);
    LinearElasticMatrix_epoxy ( 0 , 0 ) = (young_epoxy*(1.0-poisson_epoxy)/((1.0+poisson_epoxy)*(1.0-2.0*poisson_epoxy)));
    LinearElasticMatrix_epoxy ( 1 , 1 ) = LinearElasticMatrix_epoxy ( 0 , 0 );
    LinearElasticMatrix_epoxy ( 2 , 2 ) = LinearElasticMatrix_epoxy ( 0 , 0 );

    LinearElasticMatrix_epoxy ( 3 , 3 ) = LinearElasticMatrix_epoxy ( 0 , 0 )*(1.0-2.0*poisson_epoxy)/(2.0*(1.0-poisson_epoxy));
    LinearElasticMatrix_epoxy ( 4 , 4 ) = LinearElasticMatrix_epoxy ( 3 , 3 );
    LinearElasticMatrix_epoxy ( 5 , 5 ) = LinearElasticMatrix_epoxy ( 3 , 3 );

    LinearElasticMatrix_epoxy ( 0 , 1 ) = LinearElasticMatrix_epoxy ( 0 , 0 )*poisson_epoxy/(1.0-poisson_epoxy);
    LinearElasticMatrix_epoxy ( 1 , 0 ) = LinearElasticMatrix_epoxy ( 0 , 1 );

    LinearElasticMatrix_epoxy ( 0 , 2 ) = LinearElasticMatrix_epoxy ( 0 , 1 );
    LinearElasticMatrix_epoxy ( 2 , 0 ) = LinearElasticMatrix_epoxy ( 0 , 1 );

    LinearElasticMatrix_epoxy ( 1 , 2 ) = LinearElasticMatrix_epoxy ( 0 , 1 );
    LinearElasticMatrix_epoxy ( 2 , 1 ) = LinearElasticMatrix_epoxy ( 0 , 1 );

    noalias(rLinearElasticMatrix) = 0.3*LinearElasticMatrix_epoxy + 0.7*LinearElasticMatrix_glass;

}

} // Namespace Kratos
