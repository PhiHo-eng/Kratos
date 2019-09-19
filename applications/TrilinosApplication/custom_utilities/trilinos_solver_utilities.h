//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//

#if !defined (KRATOS_TRILINOS_SOLVER_UTILITIES_H_INCLUDED)
#define KRATOS_TRILINOS_SOLVER_UTILITIES_H_INCLUDED

// External includes
#include "Teuchos_ParameterList.hpp"

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"

namespace Kratos {
namespace TrilinosSolverUtilities {

void SetTeuchosParameters(const Parameters rSettings, Teuchos::ParameterList& rParameterlist);

}  // namespace TrilinosSolverUtilities.
}  // namespace Kratos.

#endif // KRATOS_TRILINOS_SOLVER_UTILITIES_H_INCLUDED defined
