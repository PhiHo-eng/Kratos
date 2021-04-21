// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//
//

/* System includes */

/* External includes */

/* Project includes */
#include "spaces/ublas_space.h"
#include "custom_strategies/custom_schemes/residual_based_relaxation_scheme.hpp"

namespace Kratos
{

///@name Type Definitions
///@{

typedef TUblasSparseSpace<double> SparseSpaceType;
typedef TUblasDenseSpace<double> LocalSpaceType;

typedef Scheme<SparseSpaceType, LocalSpaceType> SchemeType;
typedef ResidualBasedRelaxationScheme<SparseSpaceType,  LocalSpaceType> ResidualBasedRelaxationSchemeType;

//NOTE: here we must create persisting objects for the strategies
static ResidualBasedRelaxationSchemeType msResidualBasedRelaxationScheme;

template<>
std::vector<Internals::RegisteredPrototypeBase<SchemeType>> ResidualBasedRelaxationSchemeType::msPrototypes{
    Internals::RegisteredPrototype<ResidualBasedRelaxationSchemeType, SchemeType>(ResidualBasedRelaxationSchemeType::Name(), msResidualBasedRelaxationScheme)};

///@}

} /* namespace Kratos.*/