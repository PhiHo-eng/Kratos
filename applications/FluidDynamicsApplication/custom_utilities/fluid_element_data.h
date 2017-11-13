//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#if !defined(KRATOS_FLUID_ELEMENT_DATA_H)
#define KRATOS_FLUID_ELEMENT_DATA_H

// External includes

// Project includes
#include "geometries/geometry.h"
#include "includes/define.h"
#include "includes/node.h"

#include "fluid_dynamics_application_variables.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@brief Base class for data containers used within FluidElement and derived types.
template< size_t TDim, size_t TNumNodes >
class FluidElementData
{
public:
    ///@name Type Definitions
    ///@{

    using NodalScalarData = array_1d<double,TNumNodes>;

    using NodalVectorData = boost::numeric::ublas::bounded_matrix<double,TNumNodes,TDim>;

    constexpr static unsigned int Dim = TDim;

    constexpr static unsigned int NumNodes = TNumNodes;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    FluidElementData();

    /// Destructor
    virtual ~FluidElementData();

    /// (deleted) assignment operator.
    FluidElementData& operator=(FluidElementData const& rOther) = delete;

    /// (deleted) copy constructor.
    FluidElementData(FluidElementData const& rOther) = delete;

    ///@}
    ///@name Public Operations
    ///@{

    virtual void Initialize(const Element& rElement, const ProcessInfo& rProcessInfo);

    virtual void UpdateGeometryValues(double NewWeight,
        boost::numeric::ublas::matrix_row<Kratos::Matrix>& rN,
        boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim>& rDN_DX);

    ///@}
    ///@name Public Members
    ///@{
    
    double Weight;

    boost::numeric::ublas::matrix_row< Kratos::Matrix > N;

    boost::numeric::ublas::bounded_matrix< double, TNumNodes, TDim > DN_DX;
    
    ///@}
protected:

    ///@name Protected Operations
    ///@{

    void FillFromNodalData(NodalScalarData& rData, const Variable<double>& rVariable, const Geometry<Node<3>>& rGeometry);

    void FillFromNodalData(NodalVectorData& rData, const Variable<array_1d<double,3>>& rVariable, const Geometry<Node<3>>& rGeometry);

    void FillFromHistoricalNodalData(NodalScalarData& rData, const Variable<double>& rVariable, const Geometry<Node<3>>& rGeometry, const unsigned int Step);

    void FillFromHistoricalNodalData(NodalVectorData& rData, const Variable<array_1d<double,3>>& rVariable, const Geometry<Node<3>>& rGeometry, const unsigned int Step);

    void FillFromProcessInfo(double& rData, const Variable<double>& rVariable, const ProcessInfo& rProcessInfo);

    ///@}
};

///@}
}

#endif