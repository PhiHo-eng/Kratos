//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Pablo Becker
//


#if !defined(KRATOS_CONVECTION_DIFFUSION_SETTINGS_INCLUDED )
#define  KRATOS_CONVECTION_DIFFUSION_SETTINGS_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "containers/variable.h"
#include "includes/define.h"
#include "includes/kratos_components.h"


namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Convection diffusion settings. This class contains information to be used by the convection diffusion solver, all the variables that will be needed by the solver.
/** All the variables needed by any convection diffusion problem are included here. However, many problems do not require such a large ammount of variables.
 * For those cases, variables that are not defined will be set to either zero or one (depending on the varible)
 * For this purpouse, there are flags to ask wether a varible has been defined or not. For each variable, there are three main functions.
 * SetVaribleforXuse: we assing the varible for this use. When doing that we also set the flag that now this variable has been defined.
 * GetVariableforXuse: we return the variable for this use.
 * IsDefinedVariableforXuse: tells wether that variable has been defined or not.
*/
class ConvectionDiffusionSettings
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ConvectionDiffusionSettings
    KRATOS_CLASS_POINTER_DEFINITION(ConvectionDiffusionSettings);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ConvectionDiffusionSettings()
    {
        mis_defined_DensityVar=false;
        mis_defined_DiffusionVar=false;
		mis_defined_UnknownVar=false;
		mis_defined_VolumeSourceVar=false;
		mis_defined_SurfaceSourceVar=false;
		mis_defined_ProjectionVar=false;
		mis_defined_ConvectionVar=false;
		mis_defined_GradientVar = false;
		mis_defined_MeshVelocityVar=false;
		mis_defined_TransferCoefficientVar=false;
		mis_defined_VelocityVar=false;
		mis_defined_SpecificHeatVar=false;
        mis_defined_ReactionVar=false;
        mIsDefinedReactionGradientVar=false;
    };
    ConvectionDiffusionSettings(const ConvectionDiffusionSettings& rOther):
        mpDensityVar(rOther.mpDensityVar),
        mpDiffusionVar(rOther.mpDiffusionVar),
        mpUnknownVar(rOther.mpUnknownVar),
        mpVolumeSourceVar(rOther.mpVolumeSourceVar),
        mpSurfaceSourceVar(rOther.mpSurfaceSourceVar),
        mpProjectionVar(rOther. mpProjectionVar),
        mpConvectionVar(rOther.mpConvectionVar),
        mpGradientVar(rOther.mpGradientVar),
        mpMeshVelocityVar(rOther.mpMeshVelocityVar),
        mpTransferCoefficientVar(rOther.mpTransferCoefficientVar),
        mpVelocityVar(rOther.mpVelocityVar),
        mpSpecificHeatVar(rOther.mpSpecificHeatVar),
        mpReactionVar(rOther.mpReactionVar),
        mpReactionGradientVar(rOther.mpReactionGradientVar),
        mis_defined_DensityVar(rOther.mis_defined_DensityVar),
		mis_defined_DiffusionVar(rOther.mis_defined_DiffusionVar),
		mis_defined_UnknownVar(rOther.mis_defined_UnknownVar),
		mis_defined_VolumeSourceVar(rOther.mis_defined_VolumeSourceVar),
		mis_defined_SurfaceSourceVar(rOther.mis_defined_SurfaceSourceVar),
		mis_defined_ProjectionVar(rOther.mis_defined_ProjectionVar),
		mis_defined_ConvectionVar(rOther.mis_defined_ConvectionVar),
        mis_defined_GradientVar(rOther.mis_defined_GradientVar),
		mis_defined_MeshVelocityVar(rOther.mis_defined_MeshVelocityVar),
		mis_defined_TransferCoefficientVar(rOther.mis_defined_TransferCoefficientVar),
		mis_defined_VelocityVar(rOther.mis_defined_VelocityVar),
		mis_defined_SpecificHeatVar(rOther.mis_defined_SpecificHeatVar),
        mis_defined_ReactionVar(rOther.mis_defined_ReactionVar),
        mIsDefinedReactionGradientVar(rOther.mIsDefinedReactionGradientVar)
    {
    }

    /// Destructor.
    virtual ~ConvectionDiffusionSettings() {};

    ///@}
    ///@name Operators
    ///@{
    void SetDensityVariable(const Variable<double>& rvar)
    {
        mpDensityVar = &rvar;
        mis_defined_DensityVar=true;
    }
    const Variable<double>& GetDensityVariable()
    {
        return *mpDensityVar;
    }
    bool IsDefinedDensityVariable()
    {
		return mis_defined_DensityVar;
	}

    void SetDiffusionVariable(const Variable<double>& rvar)
    {
        mpDiffusionVar = &rvar;
		mis_defined_DiffusionVar=true;
    }
    const Variable<double>& GetDiffusionVariable()
    {
        return *mpDiffusionVar;
    }
    bool IsDefinedDiffusionVariable()
    {
		return mis_defined_DiffusionVar;
	}

    void SetUnknownVariable(const Variable<double>& rvar)
    {
        mpUnknownVar = &rvar;
		mis_defined_UnknownVar=true;
    }
    const Variable<double>& GetUnknownVariable()
    {
        return *mpUnknownVar;
    }
    bool IsDefinedUnknownVariable()
    {
		return mis_defined_UnknownVar;
	}

    void SetVolumeSourceVariable(const Variable<double>& rvar)
    {
        mpVolumeSourceVar = &rvar;
		mis_defined_VolumeSourceVar=true;
    }
    const Variable<double>& GetVolumeSourceVariable()
    {
        return *mpVolumeSourceVar;
    }
    bool IsDefinedVolumeSourceVariable()
    {
		return mis_defined_VolumeSourceVar;
	}

    void SetSurfaceSourceVariable(const Variable<double>& rvar)
    {
        mpSurfaceSourceVar = &rvar;
		mis_defined_SurfaceSourceVar=true;
    }
    const Variable<double>& GetSurfaceSourceVariable()
    {
        return *mpSurfaceSourceVar;
    }
    bool IsDefinedSurfaceSourceVariable()
    {
		return mis_defined_SurfaceSourceVar;
	}

    void SetProjectionVariable(const Variable<double>& rvar)
    {
        mpProjectionVar = &rvar;
		mis_defined_ProjectionVar=true;
    }
    const Variable<double>& GetProjectionVariable()
    {
        return *mpProjectionVar;
    }
    bool IsDefinedProjectionVariable()
    {
		return mis_defined_ProjectionVar;
	}

    void SetConvectionVariable(const Variable<array_1d<double,3> >& rvar)
    {
        mpConvectionVar = &rvar;
		mis_defined_ConvectionVar=true;
    }
    const Variable<array_1d<double,3> >& GetConvectionVariable()
    {
        return *mpConvectionVar;
    }
    bool IsDefinedConvectionVariable()
    {
		return mis_defined_ConvectionVar;
	}

    void SetGradientVariable(const Variable<array_1d<double,3> >& rvar)
    {
        mpGradientVar = &rvar;
        mis_defined_GradientVar=true;
    }
    const Variable<array_1d<double,3> >& GetGradientVariable() const
    {
        return *mpGradientVar;
    }
    bool IsDefinedGradientVariable() const
    {
        return mis_defined_GradientVar;
    }

    void SetMeshVelocityVariable(const Variable<array_1d<double,3> >& rvar)
    {
        mpMeshVelocityVar = &rvar;
		mis_defined_MeshVelocityVar=true;
    }
    const Variable<array_1d<double,3> >& GetMeshVelocityVariable()
    {
        return *mpMeshVelocityVar;
    }
    bool IsDefinedMeshVelocityVariable()
    {
		return mis_defined_MeshVelocityVar;
	}

    void SetTransferCoefficientVariable(const Variable<double>& rvar)
    {
        mpTransferCoefficientVar = &rvar;
		mis_defined_TransferCoefficientVar=true;
    }
    const Variable<double>& GetTransferCoefficientVariable()
    {
        return *mpTransferCoefficientVar;
    }
    bool IsDefinedTransferCoefficientVariable()
    {
		return mis_defined_TransferCoefficientVar;
	}

    void SetVelocityVariable(const Variable<array_1d<double,3> >& rvar)
    {
        mpVelocityVar = &rvar;
		mis_defined_VelocityVar=true;
    }
    const Variable<array_1d<double,3> >& GetVelocityVariable()
    {
        return *mpVelocityVar;
    }
    bool IsDefinedVelocityVariable()
    {
		return mis_defined_VelocityVar;
	}

    void SetSpecificHeatVariable(const Variable<double>& rvar)
    {
        mpSpecificHeatVar = &rvar;
		mis_defined_SpecificHeatVar=true;
    }
    const Variable<double>& GetSpecificHeatVariable()
    {
        return *mpSpecificHeatVar;
    }
    bool IsDefinedSpecificHeatVariable()
    {
		return mis_defined_SpecificHeatVar;
	}

    void SetReactionVariable(const Variable<double>& rvar)
    {
        mpReactionVar = &rvar;
		mis_defined_ReactionVar=true;
    }
    const Variable<double>& GetReactionVariable()
    {
        return *mpReactionVar;
    }
    bool IsDefinedReactionVariable()
    {
		return mis_defined_ReactionVar;
	}

    void SetReactionGradientVariable(const Variable<array_1d<double,3>>& rVar)
    {
        mpReactionGradientVar = &rVar;
		mIsDefinedReactionGradientVar=true;
    }

    const Variable<array_1d<double,3>>& GetReactionGradientVariable()
    {
        return *mpReactionGradientVar;
    }

    bool IsDefinedReactionGradientVariable()
    {
		return mIsDefinedReactionGradientVar;
	}

    ///@}
    ///@name Operations
    ///@{


    ///@}
    ///@name Access
    ///@{
    /// Assignment operator.
    ConvectionDiffusionSettings& operator=(ConvectionDiffusionSettings const& rOther)
    {
        mpDensityVar = rOther.mpDensityVar;
        mpDiffusionVar = rOther.mpDiffusionVar;
        mpUnknownVar = rOther.mpUnknownVar;
        mpVolumeSourceVar = rOther.mpVolumeSourceVar;
        mpSurfaceSourceVar = rOther.mpSurfaceSourceVar;
        mpProjectionVar = rOther.mpProjectionVar;
        mpConvectionVar = rOther.mpConvectionVar;
        mpGradientVar = rOther.mpGradientVar;
        mpMeshVelocityVar = rOther.mpMeshVelocityVar;
        mpTransferCoefficientVar = rOther.mpTransferCoefficientVar;
        mpVelocityVar = rOther.mpVelocityVar;
		mpSpecificHeatVar = rOther.mpSpecificHeatVar;
        mpReactionVar = rOther.mpReactionVar;
        mpReactionGradientVar = rOther.mpReactionGradientVar;
        //now the is_defined
        mis_defined_DensityVar = rOther.mis_defined_DensityVar;
		mis_defined_DiffusionVar = rOther.mis_defined_DiffusionVar;
		mis_defined_UnknownVar = rOther.mis_defined_UnknownVar;
		mis_defined_VolumeSourceVar = rOther.mis_defined_VolumeSourceVar;
		mis_defined_SurfaceSourceVar = rOther.mis_defined_SurfaceSourceVar;
		mis_defined_ProjectionVar = rOther.mis_defined_ProjectionVar;
        mis_defined_ConvectionVar = rOther.mis_defined_ConvectionVar;
        mis_defined_GradientVar = rOther.mis_defined_GradientVar;
		mis_defined_MeshVelocityVar = rOther.mis_defined_MeshVelocityVar;
		mis_defined_TransferCoefficientVar = rOther.mis_defined_TransferCoefficientVar;
		mis_defined_VelocityVar = rOther.mis_defined_VelocityVar;
		mis_defined_SpecificHeatVar = rOther.mis_defined_SpecificHeatVar;
        mis_defined_ReactionVar = rOther.mis_defined_ReactionVar;
        mIsDefinedReactionGradientVar = rOther.mIsDefinedReactionGradientVar;

        return *this;
    }


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "ConvectionDiffusionSettings #" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ConvectionDiffusionSettings #";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    const Variable<double>* mpDensityVar = nullptr;
    const Variable<double>* mpDiffusionVar = nullptr;
    const Variable<double>* mpUnknownVar = nullptr;
    const Variable<double>* mpVolumeSourceVar = nullptr;
    const Variable<double>* mpSurfaceSourceVar = nullptr;
    const Variable<double>* mpProjectionVar = nullptr;
    const Variable<array_1d<double,3> >* mpConvectionVar = nullptr;
    const Variable<array_1d<double,3> >* mpGradientVar = nullptr;
    const Variable<array_1d<double,3> >* mpMeshVelocityVar = nullptr;
    const Variable<double>* mpTransferCoefficientVar = nullptr;
    const Variable<array_1d<double,3> >* mpVelocityVar = nullptr;
    const Variable<double>* mpSpecificHeatVar = nullptr;
    const Variable<double>* mpReactionVar = nullptr;
    const Variable<array_1d<double,3>>* mpReactionGradientVar = nullptr;
    bool mis_defined_DensityVar = false;
    bool mis_defined_DiffusionVar = false;
    bool mis_defined_UnknownVar = false;
    bool mis_defined_VolumeSourceVar = false;
    bool mis_defined_SurfaceSourceVar = false;
    bool mis_defined_ProjectionVar = false;
    bool mis_defined_ConvectionVar = false;
    bool mis_defined_GradientVar = false;
    bool mis_defined_MeshVelocityVar = false;
    bool mis_defined_TransferCoefficientVar = false;
    bool mis_defined_VelocityVar = false;
    bool mis_defined_SpecificHeatVar = false;
    bool mis_defined_ReactionVar = false;
    bool mIsDefinedReactionGradientVar = false;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;


    virtual void save(Serializer& rSerializer) const;

    virtual void load(Serializer& rSerializer);

    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{





    ///@}

}; // Class ConvectionDiffusionSettings

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ConvectionDiffusionSettings& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ConvectionDiffusionSettings& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

#undef  KRATOS_EXPORT_MACRO
#define KRATOS_EXPORT_MACRO KRATOS_API

KRATOS_DEFINE_VARIABLE(ConvectionDiffusionSettings::Pointer, CONVECTION_DIFFUSION_SETTINGS)

#undef  KRATOS_EXPORT_MACRO
#define KRATOS_EXPORT_MACRO KRATOS_NO_EXPORT

}  // namespace Kratos.

#endif // KRATOS_CONVECTION_DIFFUSION_SETTINGS_INCLUDED  defined
