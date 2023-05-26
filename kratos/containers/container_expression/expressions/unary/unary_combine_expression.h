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

#pragma once

// System includes
#include <string>
#include <vector>

// Project includes
#include "containers/container_expression/expressions/expression.h"

namespace Kratos {

///@name Kratos Classes
///@{

class KRATOS_API(KRATOS_CORE) UnaryCombineExpression : public Expression {
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    ///@}
    ///@name Life cycle
    ///@{

    template<class TIteratorType>
    UnaryCombineExpression(
        TIteratorType Begin,
        TIteratorType End)
        : Expression(Begin != End ? (*Begin)->NumberOfEntities() : 0),
        mSourceExpressions(Begin, End)

    {
        mStrides.resize(mSourceExpressions.size());

        // every expression should have same number of entities
        IndexType local_stride = 0;
        for (IndexType i = 0; i < mSourceExpressions.size(); ++i) {
            const auto& p_expression = mSourceExpressions[i];

            KRATOS_ERROR_IF_NOT(p_expression->NumberOfEntities() == NumberOfEntities())
                << "Expression number of entities mismatch. [ required number of entities = "
                << NumberOfEntities() << ", found number of entities = "
                << p_expression->NumberOfEntities() << " ].\n"
                << "Expressions:\n"
                << "Reference = " << *mSourceExpressions[0] << "\n"
                << "Current   = " << p_expression << "\n";

            local_stride += p_expression->GetItemComponentCount();
            mStrides[i] = local_stride;
        }

        KRATOS_ERROR_IF(this->GetItemComponentCount() == 0)
            << "No expressions were given.\n";
    }

    ///@}
    ///@name Public operations
    ///@{

    template<class TIteratorType>
    static Expression::Pointer Create(
        TIteratorType Begin,
        TIteratorType End)
    {
        return Kratos::make_intrusive<UnaryCombineExpression>(Begin, End);
    }

    double Evaluate(
        const IndexType EntityIndex,
        const IndexType EntityDataBeginIndex,
        const IndexType ComponentIndex) const override;

    const std::vector<IndexType> GetItemShape() const override;

    std::string Info() const override;

    ///@}
protected:
    ///@name Private member variables
    ///@{

    const std::vector<Expression::Pointer> mSourceExpressions;

    std::vector<IndexType> mStrides;

    ///@}
};

} // namespace Kratos