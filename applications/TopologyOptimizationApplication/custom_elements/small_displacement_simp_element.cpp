//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
//                   Octaviano Malfavón Farías
//                   Eric Gonzales
//                   Philipp Hofer
//                   Erich Wehrle
//
// ==============================================================================

// Application includes

#include "small_displacement_simp_element.h"
#include "topology_optimization_application.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/constitutive_law_utilities.h"


namespace Kratos
{
SmallDisplacementSIMPElement::SmallDisplacementSIMPElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : BaseType( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

SmallDisplacementSIMPElement::SmallDisplacementSIMPElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
        : BaseType( NewId, pGeometry, pProperties )
{
    //DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer SmallDisplacementSIMPElement::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const 
{
    return Kratos::make_intrusive<SmallDisplacementSIMPElement>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer SmallDisplacementSIMPElement::Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<SmallDisplacementSIMPElement>( NewId, pGeom, pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

SmallDisplacementSIMPElement::~SmallDisplacementSIMPElement()
{
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer SmallDisplacementSIMPElement::Clone (
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    SmallDisplacementSIMPElement::Pointer p_new_elem = Kratos::make_intrusive<SmallDisplacementSIMPElement>(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    // Currently selected integration methods
    p_new_elem->SetIntegrationMethod(BaseType::mThisIntegrationMethod);

    // The vector containing the constitutive laws
    p_new_elem->SetConstitutiveLawVector(BaseType::mConstitutiveLawVector);

    return p_new_elem;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/


// =============================================================================================================================================
// STARTING / ENDING METHODS
// =============================================================================================================================================

//************************************************************************************
//************************************************************************************

	
void SmallDisplacementSIMPElement::CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Additional part for post-processing of the topology optimized model part
    if (rVariable == X_PHYS)
    {
        for (SizeType PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++)
            rValues[PointNumber] = this->GetValue(X_PHYS);
    }
    else {

            SmallDisplacement::CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
      	}

    KRATOS_CATCH( "" )
    } 

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementSIMPElement::Calculate(const Variable<double> &rVariable, double &rOutput, const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rVariable == DCDX || rVariable == DCDX_COMPLIANT || rVariable == LOCAL_STRAIN_ENERGY || rVariable == LOCAL_STRAIN_ENERGY_COMPLIANT || rVariable == DCDX_DISPLACMENT_CONTROLLED_X  || rVariable == DCDX_DISPLACMENT_CONTROLLED_Y  || rVariable == DCDX_DISPLACMENT_CONTROLLED_Z)
    {
        // Get values
        const double E_min     = this->GetProperties()[YOUNGS_MODULUS_MIN];
        const double E_initial = this->GetProperties()[YOUNGS_MODULUS_0];
        const double E_current = this->GetValue(YOUNG_MODULUS);
        const double penalty   = this->GetValue(PENAL);
        const double x_phys    = this->GetValue(X_PHYS);

        // Get element stiffness matrix and modify it with the factor that considers the adjusted Youngs Modulus according to the SIMP method
        // Note that Ke0 is computed based on the originally provided Youngs Modulus in the .mdpa-file
        MatrixType Ke0 = Matrix();
        this->CalculateLeftHandSide(Ke0, const_cast <ProcessInfo&>(rCurrentProcessInfo));
        double E_new     = (E_min + pow(x_phys, penalty) * (E_initial - E_min));

        ///Normalize the youngs modulus
        double factor    = E_new/E_current; 
        MatrixType Ke = Ke0 * factor;

        // Loop through nodes of elements and create elemental displacement vector "ue"
        Element::GeometryType& rGeom = this->GetGeometry();
        unsigned int NumNodes = rGeom.PointsNumber(); 

        // Resize "ue" according to element type
        Vector ue;
        ue.resize(NumNodes * 3);

        // Set the displacement obtained from the FE-Analysis
        for (unsigned int node_i = 0; node_i < NumNodes; node_i++) {
            array_1d<double, 3> &CurrentDisplacement = rGeom[node_i].FastGetSolutionStepValue(DISPLACEMENT);
            ue[3 * node_i + 0] = CurrentDisplacement[0];
            ue[3 * node_i + 1] = CurrentDisplacement[1];
            ue[3 * node_i + 2] = CurrentDisplacement[2];
            }

        

        if (rVariable == DCDX_COMPLIANT || rVariable == LOCAL_STRAIN_ENERGY_COMPLIANT )
        {

            Vector lambda;
            lambda = this->GetValue(LAMBDA_ADJOINT);
            //KRATOS_INFO("[TopOpt]") << "  LAMBDA - U " << lambda << " ] " << std::endl;
            //lambda = -lambda;
            Vector intermediateVector;
            intermediateVector.resize(NumNodes * 3);
            intermediateVector = prod(trans(lambda), Ke0);
            double lambda_Ke0_ue = inner_prod(intermediateVector, ue);
            if ( std::isnan(lambda_Ke0_ue))
            {
                KRATOS_ERROR << "The nan comes from the Element in TopOPt: " <<lambda <<"and "<< x_phys<<std::endl;
            }
            
            if (rVariable == DCDX_COMPLIANT)
            {
                // Calculation of the compliance sensitivities DCDX
                double dcdx = (penalty)* (E_initial - E_min) * pow(x_phys, penalty - 1) * lambda_Ke0_ue;
                this->SetValue(DCDX_COMPLIANT, dcdx); 
                //KRATOS_INFO("[TopOpt]") << "  SENSITIVITY " << ue << " ] " << std::endl;
            }
            if (rVariable == LOCAL_STRAIN_ENERGY_COMPLIANT)
            {
                // Calculation of the local strain energy (requires Ke)
                double local_strain_energy =  factor * lambda_Ke0_ue;
                this->SetValue(LOCAL_STRAIN_ENERGY_COMPLIANT, local_strain_energy);

            }

        }

        else if (rVariable == DCDX || rVariable == LOCAL_STRAIN_ENERGY)
        {
            

            // Calculate trans(ue)*Ke0*ue
            Vector intermediateVector;
            intermediateVector.resize(NumNodes * 3);
            intermediateVector = prod(trans(ue), Ke0);
            double ue_Ke0_ue = inner_prod(intermediateVector, ue);

            if (rVariable == DCDX)
            {
                // Calculation of the compliance sensitivities DCDX
                double dcdx = (-penalty)* (E_initial - E_min) * pow(x_phys, penalty - 1) * ue_Ke0_ue;
                this->SetValue(DCDX, dcdx); 
            }
            if (rVariable == LOCAL_STRAIN_ENERGY)
            {
                // Calculation of the local strain energy (requires Ke)
                double local_strain_energy = factor * ue_Ke0_ue;
                this->SetValue(LOCAL_STRAIN_ENERGY, local_strain_energy);

            }
        }

        else if (rVariable == DCDX_DISPLACMENT_CONTROLLED_X )
        {
            Vector lambda;
            Vector lambda_controlled;
            lambda = this->GetValue(LAMBDA_ADJOINT);

            for (unsigned int node_i = 0; node_i < NumNodes; node_i++) {
            array_1d<double, 3> &CurrentDisplacement = rGeom[node_i].FastGetSolutionStepValue(DISPLACEMENT);
            ue[3 * node_i + 0] = CurrentDisplacement[0];
            lambda_controlled[3* node_i +0] = lambda[3*node_i+0];
            }

            Vector intermediateVector;
            intermediateVector.resize(NumNodes * 3);
            intermediateVector = prod(trans(ue), Ke0);
            double ue_Ke0_ue = inner_prod(intermediateVector, ue);

            if (rVariable == DCDX_DISPLACMENT_CONTROLLED_X)
            {
                // Calculation of the compliance sensitivities DCDX
                double dcdx = (-penalty)* (E_initial - E_min) * pow(x_phys, penalty - 1) * ue_Ke0_ue;
                this->SetValue(DCDX_DISPLACMENT_CONTROLLED_X, dcdx);
            }
        }

        else if (rVariable == DCDX_DISPLACMENT_CONTROLLED_Y )
        {
            Vector lambda;
            Vector lambda_controlled;
            lambda = this->GetValue(LAMBDA_ADJOINT);

            for (unsigned int node_i = 0; node_i < NumNodes; node_i++) {
            array_1d<double, 3> &CurrentDisplacement = rGeom[node_i].FastGetSolutionStepValue(DISPLACEMENT);
            ue[3 * node_i + 1] = CurrentDisplacement[1];
            lambda_controlled[3* node_i + 1] = lambda[3*node_i+1];
            }

            Vector intermediateVector;
            intermediateVector.resize(NumNodes * 3);
            intermediateVector = prod(trans(ue), Ke0);
            double ue_Ke0_ue = inner_prod(intermediateVector, ue);

            if (rVariable == DCDX_DISPLACMENT_CONTROLLED_Y)
            {
                // Calculation of the compliance sensitivities DCDX
                double dcdx = (-penalty)* (E_initial - E_min) * pow(x_phys, penalty - 1) * ue_Ke0_ue;
                this->SetValue(DCDX_DISPLACMENT_CONTROLLED_Y, dcdx);
            }
        }

        else if (rVariable == DCDX_DISPLACMENT_CONTROLLED_Z )
        {
            Vector lambda;
            Vector lambda_controlled;
            lambda = this->GetValue(LAMBDA_ADJOINT);

            for (unsigned int node_i = 0; node_i < NumNodes; node_i++) {
            array_1d<double, 3> &CurrentDisplacement = rGeom[node_i].FastGetSolutionStepValue(DISPLACEMENT);
            ue[3 * node_i + 0] = CurrentDisplacement[2];
            lambda_controlled[3* node_i +2] = lambda[3*node_i+2];
            }

            Vector intermediateVector;
            intermediateVector.resize(NumNodes * 3);
            intermediateVector = prod(trans(ue), Ke0);
            double ue_Ke0_ue = inner_prod(intermediateVector, ue);

            if (rVariable == DCDX_DISPLACMENT_CONTROLLED_Z)
            {
                // Calculation of the compliance sensitivities DCDX
                double dcdx = (-penalty)* (E_initial - E_min) * pow(x_phys, penalty - 1) * ue_Ke0_ue;
                this->SetValue(DCDX_DISPLACMENT_CONTROLLED_Z, dcdx);
            }
        }




    } else if (rVariable == DVDX) {
		// Calculation of the volume sensitivities DVDX
        double element_size = this->GetValue(INITIAL_ELEMENT_SIZE);
        if ( std::isnan(element_size))
        if ( std::isnan(element_size))
        {
            KRATOS_ERROR << "The nan comes from the Element in TopOpt from element size: " <<element_size<<std::endl;
        }
        this->SetValue(DVDX, element_size*1);
	}

    KRATOS_CATCH( "" )
} 
 
/***********************************************************************************/
/***********************************************************************************/


void SmallDisplacementSIMPElement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType )
}

void SmallDisplacementSIMPElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType )
}

} // Namespace Kratos