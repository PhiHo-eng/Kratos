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

#if !defined (KRATOS_JOINT_BILINEAR_COHESIVE_3D_LAW_H_INCLUDED)
#define  KRATOS_JOINT_BILINEAR_COHESIVE_3D_LAW_H_INCLUDED

// System includes
#include <cmath>

// Project includes
#include "includes/serializer.h"
#include "includes/checks.h"
#include "includes/constitutive_law.h"

// Application includes
#include "dam_application_variables.h"
#include "poromechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(DAM_APPLICATION) JointBilinearCohesive3DLaw : public ConstitutiveLaw
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(JointBilinearCohesive3DLaw);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default Constructor
    JointBilinearCohesive3DLaw()
    {
    }

    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<JointBilinearCohesive3DLaw>(JointBilinearCohesive3DLaw(*this));
    }

    // Copy Constructor
    JointBilinearCohesive3DLaw (const JointBilinearCohesive3DLaw& rOther) : ConstitutiveLaw(rOther)
    {
    }

    // Destructor
    ~JointBilinearCohesive3DLaw() override
    {
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void GetLawFeatures(Features& rFeatures) override;

    int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo) const override;

    void InitializeMaterial( const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const Vector& rShapeFunctionsValues ) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateMaterialResponseCauchy (Parameters & rValues) override;

    void FinalizeMaterialResponseCauchy (Parameters & rValues) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    double& GetValue( const Variable<double>& rThisVariable, double& rValue ) override;

    void SetValue( const Variable<double>& rVariable, const double& rValue, const ProcessInfo& rCurrentProcessInfo ) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    struct ConstitutiveLawVariables
    {
        double CriticalDisplacement;
        double CriticalDisplacementTangent;
        double DamageThreshold;
        double YieldStress;
        double YoungModulus;
        double PoissonCoefficient;
        double FrictionCoefficient;
        double PenaltyStiffness;
        double MaxTensileStress;
        double MaxCompresiveStress;
        double Cohesion;
        double YieldStressTangent;

        Matrix CompressionMatrix;
        Matrix WeightMatrix;

        double EquivalentStrain;
        bool LoadingFlag;
        double LoadingFunction;
    };

    // Member Variables

    double mStateVariable;
    double mUpliftPressure;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    virtual void InitializeConstitutiveLawVariables(ConstitutiveLawVariables& rVariables, Parameters& rValues);

    virtual void ComputeEquivalentStrain(ConstitutiveLawVariables& rVariables, Parameters& rValues);

    virtual void CheckLoadingFunction(ConstitutiveLawVariables& rVariables, Parameters& rValues);

    virtual void ComputeConstitutiveMatrix(Matrix& rConstitutiveMatrix,
                                            ConstitutiveLawVariables& rVariables,
                                            Parameters& rValues);

    virtual void ComputeStressVector(Vector& rStressVector,
                                        ConstitutiveLawVariables& rVariables,
                                        Parameters& rValues);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    // Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw )
    }

}; // Class JointBilinearCohesive3DLaw
}  // namespace Kratos.
#endif // KRATOS_JOINT_BILINEAR_COHESIVE_3D_LAW_H_INCLUDED  defined
