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
#include <variant>
#include <string>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/kratos_parameters.h"
#include "controllers/controller.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/**
 * @class Temporal controller
 * @ingroup KratosCore
 * @brief A temporal controller to control behavior based on temporal values.
 * @details This controller checks prescribed control variables, and evalutes.
 * @author Suneth Warnakulasuriya
*/
class TemporalController: public Controller
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TemporalController
    KRATOS_CLASS_POINTER_DEFINITION(TemporalController);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    TemporalController(
        const Model& rModel,
        Parameters Settings);

    /// Destructor.
    ~TemporalController() override = default;

    /// Copy constructor.
    //TODO: Check. It is required by the registry
    TemporalController(TemporalController const& rOther) = default;

    /// Move constructor
    TemporalController(TemporalController&& rOther) noexcept = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    TemporalController& operator=(TemporalController const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method creates an pointer of the controller
     * @details We consider as input, a Model and a set of Parameters for the sake of generality
     * @warning Must be overrided in controller implementation
     * @param rModel The model to be considered
     * @param ThisParameters The configuration parameters
     */
    Controller::Pointer Create(
        Model& rModel,
        Parameters Settings) const override;

    /**
     * @brief Checks that input conditions are correct.
     */
    int Check() const override;

    /**
     * @brief Using input data, returns bool.
     */
    bool Evaluate() override;

    /**
     * @brief Get the current control value.
    */
    std::variant<int, double> GetCurrentControlValue() const;

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    Parameters GetDefaultParameters() const override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

    ///@}

private:
    ///@name Private member variables
    ///@{

    Model const * mpModel;

    std::variant<int, double> mNextOutput;

    std::variant<int, double> mInterval;

    std::variant<
        const Variable<int>*,
        const Variable<double>*> mpVariable;

    std::string mModelPartName;

    ///@}
    ///@name Private operations
    ///@{

    void ScheduleNextOutput();

    ///@}

}; // Class TemporalController

}  // namespace Kratos.
