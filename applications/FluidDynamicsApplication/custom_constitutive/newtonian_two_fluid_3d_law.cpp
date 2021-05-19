//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Daniel Diez
//

// System includes

// External includes

// Project includes
#include "includes/cfd_variables.h"
#include "includes/checks.h"
#include "utilities/element_size_calculator.h"

// Application includes
#include "custom_constitutive/newtonian_two_fluid_3d_law.h"
#include "fluid_dynamics_application_variables.h"
namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

NewtonianTwoFluid3DLaw::NewtonianTwoFluid3DLaw()
    : Newtonian3DLaw()
{}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

NewtonianTwoFluid3DLaw::NewtonianTwoFluid3DLaw(const NewtonianTwoFluid3DLaw& rOther)
    : Newtonian3DLaw(rOther)
{}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer NewtonianTwoFluid3DLaw::Clone() const {
    return Kratos::make_shared<NewtonianTwoFluid3DLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

NewtonianTwoFluid3DLaw::~NewtonianTwoFluid3DLaw() {}


std::string NewtonianTwoFluid3DLaw::Info() const {
    return "NewtonianTwoFluid3DLaw";
}

int NewtonianTwoFluid3DLaw::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo)
{
    for (unsigned int i = 0; i < rElementGeometry.size(); i++) {
        const Node<3>& rNode = rElementGeometry[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DYNAMIC_VISCOSITY,rNode);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DENSITY,rNode);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISTANCE,rNode);
        KRATOS_ERROR_IF(rNode.GetSolutionStepValue(DYNAMIC_VISCOSITY) <= 0.0)
            << "DYNAMIC_VISCOSITY was not correctly assigned to the nodes for Constitutive Law.\n";
        KRATOS_ERROR_IF(rNode.GetSolutionStepValue(DENSITY) <= 0.0)
            << "DENSITY was not correctly assigned to the nodes for Constitutive Law.\n";
    }
    return 0;
}


double NewtonianTwoFluid3DLaw::GetEffectiveViscosity(ConstitutiveLaw::Parameters& rParameters) const
{
    double viscosity;
    double artificial_viscosity;
    EvaluateInPoint(viscosity, DYNAMIC_VISCOSITY, rParameters);
    const Properties& prop = rParameters.GetMaterialProperties();
    EvaluateInPointShockCapturing(artificial_viscosity,rParameters);
    if (prop.Has(C_SMAGORINSKY)) {
        const double csmag = prop[C_SMAGORINSKY];
        if (csmag > 0.0) {
            double density;
            EvaluateInPoint(density, DENSITY, rParameters);
            const double strain_rate = EquivalentStrainRate(rParameters);
            const BoundedMatrix<double, 4, 3>& DN_DX = rParameters.GetShapeFunctionsDerivatives();
            const double elem_size = ElementSizeCalculator<3,4>::GradientsElementSize(DN_DX);
            const double length_scale = std::pow(csmag * elem_size, 2);
            // length_scale *= length_scale;
            viscosity += 2.0*length_scale * strain_rate * density;
        }
    }
    // Implementing shock capturing option for two fluid in 3d law
    // Ask RUben 
  
        
    double added_shock_capturing_viscosity= artificial_viscosity+viscosity;

    return added_shock_capturing_viscosity;
}


void NewtonianTwoFluid3DLaw::EvaluateInPoint(double& rResult,
    const Variable<double>& rVariable,
    ConstitutiveLaw::Parameters& rParameters) const {

    const SizeType nnodes = 4;
    const GeometryType& geom = rParameters.GetElementGeometry();
    const array_1d<double,nnodes>& N = rParameters.GetShapeFunctionsValues();

    //compute sign of distance on gauss point
    double dist = 0.0;
    for (unsigned int i = 0; i < nnodes; i++)
        dist += N[i] * geom[i].FastGetSolutionStepValue(DISTANCE);

    SizeType navg = 0; //number of nodes on the same side as the gauss point
    double value = 0.0;
    for (unsigned int i = 0; i < nnodes; i++) {
        if ( dist * geom[i].FastGetSolutionStepValue(DISTANCE) > 0.0) {
            navg += 1;
            value += geom[i].FastGetSolutionStepValue(rVariable);
        }
    }
    rResult = value/navg;
}



void NewtonianTwoFluid3DLaw::EvaluateInPointShockCapturing(double& rResult,
    ConstitutiveLaw::Parameters& rParameters) const {

    const SizeType nnodes = 4;
    const GeometryType& geom = rParameters.GetElementGeometry();
    const array_1d<double,nnodes>& N = rParameters.GetShapeFunctionsValues();

    //compute sign of distance on gauss point
    
    double aritficial_viscosity_shockcapturing=0;
    SizeType navg = 0; //number of nodes on the same side as the gauss point
    
    for (unsigned int i = 0; i < nnodes; i++){
        
    
        navg += 1;
        aritficial_viscosity_shockcapturing += N[i] * geom[i].GetValue(ARTIFICIAL_DYNAMIC_VISCOSITY);
    }
    rResult = aritficial_viscosity_shockcapturing/navg;
    }
 
    
    




double NewtonianTwoFluid3DLaw::EquivalentStrainRate(ConstitutiveLaw::Parameters& rParameters) const {

    const Vector& S = rParameters.GetStrainVector();

    // Norm of symetric gradient (cross terms don't get the 2)
    return std::sqrt(2.*S[0]*S[0] + 2.*S[1]*S[1] + 2.*S[2]*S[2] + S[3]*S[3] + S[4]*S[4] + S[5]*S[5]);
}



void NewtonianTwoFluid3DLaw::save(Serializer& rSerializer) const {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Newtonian3DLaw )
}

void NewtonianTwoFluid3DLaw::load(Serializer& rSerializer) {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Newtonian3DLaw )
}

} // Namespace Kratos
