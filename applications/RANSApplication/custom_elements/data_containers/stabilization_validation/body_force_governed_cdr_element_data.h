//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#if !defined(KRATOS_STABILIZATION_VALIDATION_BODY_FORCE_GOVERNED_CDR_ELEMENT_DATA_H)
#define KRATOS_STABILIZATION_VALIDATION_BODY_FORCE_GOVERNED_CDR_ELEMENT_DATA_H

// System includes

// Project includes
#include "containers/variable.h"
#include "geometries/geometry_data.h"
#include "includes/node.h"
#include "includes/process_info.h"
#include "includes/ublas_interface.h"

// Application includes
#include "custom_elements/convection_diffusion_reaction_element_data.h"

namespace Kratos
{
///@name  Functions
///@{

namespace StabilizationValidation
{
class BodyForceGovernedCDRElementData : public ConvectionDiffusionReactionElementData<2>
{
public:
    ///@name Type Definitions
    ///@{

    using BaseType = ConvectionDiffusionReactionElementData<2>;

    using GeometryType = typename BaseType::GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    BodyForceGovernedCDRElementData(
        const GeometryType& rGeometry,
        const Properties& rProperties,
        const ProcessInfo& rProcessInfo);

    ~BodyForceGovernedCDRElementData() override = default;

    ///@}
    ///@name Static Operations
    ///@{

    static const Variable<double>& GetScalarVariable();

    static void Check(
        const Element& rElement,
        const ProcessInfo& rCurrentProcessInfo);

    static const std::string GetName()
    {
        return "BodyForceGovernedCDRElementData";
    }

    ///@}
    ///@name Operations

    void CalculateConstants(
        const ProcessInfo& rCurrentProcessInfo);

    void CalculateGaussPointData(
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives,
        const int Step = 0);

    ///@}

protected:
    ///@name Protected Members
    ///@{

    using BaseType::mEffectiveVelocity;
    using BaseType::mEffectiveKinematicViscosity;
    using BaseType::mReactionTerm;
    using BaseType::mSourceTerm;

    BoundedMatrix<double, 2, 3> mNodalCoordinates;

    ///@}
};

///@}

} // namespace StabilizationValidation

} // namespace Kratos

#endif // KRATOS_STABILIZATION_VALIDATION_BODY_FORCE_GOVERNED_CDR_ELEMENT_DATA_H