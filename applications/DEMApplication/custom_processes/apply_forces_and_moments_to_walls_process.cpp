#include "apply_forces_and_moments_to_walls_process.h"

namespace Kratos
{
    /* Public functions *******************************************************/
    ApplyForcesAndMomentsToWallsProcess::ApplyForcesAndMomentsToWallsProcess(
        ModelPart& rModelPart,
        Parameters rParameters
        ) : Process(Flags()) , mrModelPart(rModelPart), mParameters(rParameters), mInterval(rParameters)
    {
        KRATOS_TRY

        //only include validation with c++11 since raw_literals do not exist in c++03
        Parameters default_parameters( R"(
            {
                "help"                 : "This process applies constraints to the particles in a certain submodelpart, for a certain time interval",
                "mesh_id"              : 0,
                "model_part_name"      : "please_specify_model_part_name",
                "force_settings" : {
                    "value"                : [10.0, "3*t", "x+y"]
                },
                "moment_settings" : {
                    "value"                : [10.0, "3*t", "x+y"]
                },
                "interval"             : [0.0, 1e30]
            } )" );



        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["model_part_name"];

        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mForceFunctions.clear();
        mMomentFunctions.clear();

        for(int i=0; i<3; i++) {
            if(rParameters["force_settings"]["value"][i].IsNumber()) {
                mForceValueIsNumeric[i] = true;
                mForceValues[i] = rParameters["force_settings"]["value"][i].GetDouble();
                mForceFunctions.push_back(PythonGenericFunctionUtility("0.0")); // because I can't construct an array_1d of these
            }
            else {
                mForceValueIsNumeric[i] = false;
                mForceFunctions.push_back(PythonGenericFunctionUtility(rParameters["force_settings"]["value"][i].GetString()));
            }
            if(rParameters["moment_settings"]["value"][i].IsNumber()) {
                mMomentValueIsNumeric[i] = true;
                mMomentValues[i] = rParameters["moment_settings"]["value"][i].GetDouble();
                mMomentFunctions.push_back(PythonGenericFunctionUtility("0.0")); // because I can't construct an array_1d of these
            }
            else {
                mMomentValueIsNumeric[i] = false;
                mMomentFunctions.push_back(PythonGenericFunctionUtility(rParameters["moment_settings"]["value"][i].GetString()));
            }
        }

        mParameters = rParameters;

        KRATOS_CATCH("");
    }

    ApplyForcesAndMomentsToWallsProcess::~ApplyForcesAndMomentsToWallsProcess()
    {

    }

    void ApplyForcesAndMomentsToWallsProcess::Execute()
    {
    }

    void ApplyForcesAndMomentsToWallsProcess::ExecuteInitializeSolutionStep()
    {
        KRATOS_TRY;


        const double time = mrModelPart.GetProcessInfo()[TIME];

        if(!mInterval.IsInInterval(time)) return;

        for(int i=0; i<3; i++) {
            double force_value = 0.0;
            if(mForceValueIsNumeric[i]) {
                force_value = mForceValues[i];
            }
            else {
                force_value = mForceFunctions[i].CallFunction(mrModelPart[RIGID_BODY_CENTER_OF_MASS][0],
                                                              mrModelPart[RIGID_BODY_CENTER_OF_MASS][1],
                                                              mrModelPart[RIGID_BODY_CENTER_OF_MASS][2],
                                                              time);
            }
            mrModelPart[EXTERNAL_APPLIED_FORCE][i] = force_value;

            double moment_value = 0.0;
            if(mMomentValueIsNumeric[i]) {
                moment_value = mMomentValues[i];
            }
            else {
                    moment_value = mMomentFunctions[i].CallFunction(mrModelPart[RIGID_BODY_CENTER_OF_MASS][0],
                                                                    mrModelPart[RIGID_BODY_CENTER_OF_MASS][1],
                                                                    mrModelPart[RIGID_BODY_CENTER_OF_MASS][2],
                                                                    time);
            }
            mrModelPart[EXTERNAL_APPLIED_MOMENT][i] = moment_value;
        }

        KRATOS_CATCH("");
    }

    void ApplyForcesAndMomentsToWallsProcess::ExecuteFinalizeSolutionStep()
    {
        KRATOS_TRY;


        const double time = mrModelPart.GetProcessInfo()[TIME];

        if(mInterval.IsInInterval(time)) return;

        mrModelPart[EXTERNAL_APPLIED_FORCE] = ZeroVector(3);
        mrModelPart[EXTERNAL_APPLIED_MOMENT] = ZeroVector(3);

        KRATOS_CATCH("");
    }

    std::string ApplyForcesAndMomentsToWallsProcess::Info() const
    {
        return "ApplyForcesAndMomentsToWallsProcess";
    }

    void ApplyForcesAndMomentsToWallsProcess::PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ApplyForcesAndMomentsToWallsProcess";
    }

    void ApplyForcesAndMomentsToWallsProcess::PrintData(std::ostream& rOStream) const
    {
    }

    /* External functions *****************************************************/

    /// output stream function
    inline std::ostream& operator << (
        std::ostream& rOStream,
        const ApplyForcesAndMomentsToWallsProcess& rThis)
    {
        rThis.PrintData(rOStream);
        return rOStream;
    }

}
