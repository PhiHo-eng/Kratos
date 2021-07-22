//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

#if !defined(KRATOS_INTEGRATION_INFO_H_INCLUDED )
#define  KRATOS_INTEGRATION_INFO_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "geometries/geometry_data.h"
#include "containers/data_value_container.h"
#include "containers/flags.h"

#include "integration_flags.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/// Integration information for the creation of integration points.
/* Within this class distinct information of integration can be
 * stored and processed.
 */
class IntegrationInfo : public Flags
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of IntegrationPoint
    KRATOS_CLASS_POINTER_DEFINITION(IntegrationInfo);

    typedef typename Point::IndexType SizeType;
    typedef typename Point::IndexType IndexType;

    /// Integration methods implemented specified within enum.
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    ///@}
    ///@name Local Flags
    ///@{

    KRATOS_DEFINE_LOCAL_FLAG(DO_NOT_CREATE_TESSELLATION_ON_SLAVE);

    ///@}
    ///@name Type Definitions
    ///@{

    enum class QuadratureMethod
    {
        Default,
        GAUSS,
        EXTENDED_GAUSS
    };

    ///@}
    ///@name Life Cycle
    ///@{

    IntegrationInfo(SizeType LocalSpaceDimension,
        IntegrationMethod ThisIntegrationMethod);

    IntegrationInfo(SizeType LocalSpaceDimension,
        SizeType NumberOfIntegrationPointsPerSpan,
        QuadratureMethod ThisQuadratureMethod = QuadratureMethod::GAUSS);

    IntegrationInfo(
        std::vector<SizeType> NumberOfIntegrationPointsPerSpanVector,
        std::vector<QuadratureMethod> ThisQuadratureMethodVector);

    ///@}
    ///@name Dimension
    ///@{

    SizeType LocalSpaceDimension()
    {
        return mNumberOfIntegrationPointsPerSpanVector.size();
    }

    ///@}
    ///@name integration rules
    ///@{

    void SetIntegrationMethod(
        IndexType DimensionIndex,
        IntegrationMethod ThisIntegrationMethod);

    SizeType GetNumberOfIntegrationPointsPerSpan(IndexType DimensionIndex) const;

    void SetNumberOfIntegrationPointsPerSpan(IndexType DimensionIndex,
        SizeType NumberOfIntegrationPointsPerSpan);

    QuadratureMethod GetQuadratureMethodVector(IndexType DimensionIndex) const;

    void SetQuadratureMethodVector(IndexType DimensionIndex,
        QuadratureMethod ThisQuadratureMethod);

    /* returns the IntegrationMethod to
     * corresponding to the direction index.
     */
    IntegrationMethod GetIntegrationMethod(
        IndexType DimensionIndex) const;

    /* Evaluates the corresponding IntegrationMethod to 
     * the number of points and the quadrature method.
     */
    static IntegrationMethod GetIntegrationMethod(
        SizeType NumberOfIntegrationPointsPerSpan,
        QuadratureMethod ThisQuadratureMethod);

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << " Integration info with local space dimension: " << mNumberOfIntegrationPointsPerSpanVector.size()
            << " and number of integration points per spans: " << mNumberOfIntegrationPointsPerSpanVector;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << " Integration info with local space dimension: " << mNumberOfIntegrationPointsPerSpanVector.size()
            << " and number of integration points per spans: " << mNumberOfIntegrationPointsPerSpanVector;
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << " Integration info with local space dimension: " << mNumberOfIntegrationPointsPerSpanVector.size()
            << " and number of integration points per spans: " << mNumberOfIntegrationPointsPerSpanVector;
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    std::vector<SizeType> mNumberOfIntegrationPointsPerSpanVector;

    std::vector<QuadratureMethod> mQuadratureMethodVector;

    ///@}

}; // Class IntegrationPoint

///@}
///@name Input and output
///@{


/// input stream function
template<std::size_t TDimension, class TDataType, class TWeightType>
inline std::istream& operator >> (std::istream& rIStream, IntegrationInfo& rThis);

/// output stream function
template<std::size_t TDimension, class TDataType, class TWeightType>
inline std::ostream& operator << (std::ostream& rOStream, const IntegrationInfo& rThis)
{
    rThis.PrintInfo(rOStream);
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_INTEGRATION_INFO_H_INCLUDED  defined