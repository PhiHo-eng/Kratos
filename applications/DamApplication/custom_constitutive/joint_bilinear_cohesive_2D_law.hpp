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

#if !defined (KRATOS_JOINT_BILINEAR_COHESIVE_2D_LAW_H_INCLUDED)
#define  KRATOS_JOINT_BILINEAR_COHESIVE_2D_LAW_H_INCLUDED

// Application includes
#include "custom_constitutive/joint_bilinear_cohesive_3D_law.hpp"

namespace Kratos
{

class KRATOS_API(DAM_APPLICATION) JointBilinearCohesive2DLaw : public JointBilinearCohesive3DLaw
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(JointBilinearCohesive2DLaw);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default Constructor
    JointBilinearCohesive2DLaw()
    {
    }

    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<JointBilinearCohesive2DLaw>(JointBilinearCohesive2DLaw(*this));
    }

    // Copy Constructor
    JointBilinearCohesive2DLaw (const JointBilinearCohesive2DLaw& rOther) : JointBilinearCohesive3DLaw(rOther)
    {
    }

    // Destructor
    ~JointBilinearCohesive2DLaw() override
    {
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void GetLawFeatures(Features& rFeatures) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ComputeEquivalentStrain(ConstitutiveLawVariables& rVariables, Parameters& rValues) override;

    void ComputeConstitutiveMatrix(Matrix& rConstitutiveMatrix,
                                    ConstitutiveLawVariables& rVariables,
                                    Parameters& rValues) override;

    void ComputeStressVector(Vector& rStressVector,
                                ConstitutiveLawVariables& rVariables,
                                Parameters& rValues) override;

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

}; // Class JointBilinearCohesive2DLaw
}  // namespace Kratos.
#endif // KRATOS_JOINT_BILINEAR_COHESIVE_2D_LAW_H_INCLUDED  defined
