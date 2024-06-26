//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:              LMonforte $
//   Date:                $Date:                July 2015 $
//   Revision:            $Revision:                 -0.1 $
//
//

// System includes

// External includes

// Project includes
#include "custom_elements/updated_lagrangian_U_J_wP_element.hpp"
#include "pfem_solid_mechanics_application_variables.h"

#include "custom_utilities/water_pressure_utilities_Jacobian.hpp"

namespace Kratos
{


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************
   // Aquest a l'altre no hi és....
   UpdatedLagrangianUJwPElement::UpdatedLagrangianUJwPElement()
      : UpdatedLagrangianUJElement()
   {
      //DO NOT CALL IT: only needed for Register and Serialization!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   UpdatedLagrangianUJwPElement::UpdatedLagrangianUJwPElement( IndexType NewId, GeometryType::Pointer pGeometry )
      : UpdatedLagrangianUJElement( NewId, pGeometry )
   {
      //DO NOT ADD DOFS HERE!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   UpdatedLagrangianUJwPElement::UpdatedLagrangianUJwPElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
      : UpdatedLagrangianUJElement( NewId, pGeometry, pProperties )
   {
   }


   //******************************COPY CONSTRUCTOR**************************************
   //************************************************************************************

   UpdatedLagrangianUJwPElement::UpdatedLagrangianUJwPElement( UpdatedLagrangianUJwPElement const& rOther)
      :UpdatedLagrangianUJElement(rOther)
       , mTimeStep(rOther.mTimeStep)
   {
   }


   //*******************************ASSIGMENT OPERATOR***********************************
   //************************************************************************************

   UpdatedLagrangianUJwPElement&  UpdatedLagrangianUJwPElement::operator=(UpdatedLagrangianUJwPElement const& rOther)
   {
      UpdatedLagrangianUJElement::operator=(rOther);

      mTimeStep = rOther.mTimeStep;

      return *this;
   }


   //*********************************OPERATIONS*****************************************
   //************************************************************************************

   Element::Pointer UpdatedLagrangianUJwPElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
   {
      return Element::Pointer( new UpdatedLagrangianUJwPElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
   }


   //************************************CLONE*******************************************
   //************************************************************************************

   Element::Pointer UpdatedLagrangianUJwPElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
   {

      UpdatedLagrangianUJwPElement NewElement( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

      //-----------//

      NewElement.mThisIntegrationMethod = mThisIntegrationMethod;


      if ( NewElement.mConstitutiveLawVector.size() != mConstitutiveLawVector.size() )
      {
         NewElement.mConstitutiveLawVector.resize(mConstitutiveLawVector.size());

         if( NewElement.mConstitutiveLawVector.size() != NewElement.GetGeometry().IntegrationPointsNumber() )
            KRATOS_THROW_ERROR( std::logic_error, "constitutive law not has the correct size ", NewElement.mConstitutiveLawVector.size() )
      }

      for(unsigned int i=0; i<mConstitutiveLawVector.size(); i++)
      {
         NewElement.mConstitutiveLawVector[i] = mConstitutiveLawVector[i]->Clone();
      }

      //-----------//

      if ( NewElement.mDeformationGradientF0.size() != mDeformationGradientF0.size() )
         NewElement.mDeformationGradientF0.resize(mDeformationGradientF0.size());

      for(unsigned int i=0; i<mDeformationGradientF0.size(); i++)
      {
         NewElement.mDeformationGradientF0[i] = mDeformationGradientF0[i];
      }

      NewElement.mDeterminantF0 = mDeterminantF0;

      NewElement.mTimeStep = mTimeStep;

      NewElement.SetData(this->GetData());
      NewElement.SetFlags(this->GetFlags());

      return Element::Pointer( new UpdatedLagrangianUJwPElement(NewElement) );
   }


   //*******************************DESTRUCTOR*******************************************
   //************************************************************************************

   UpdatedLagrangianUJwPElement::~UpdatedLagrangianUJwPElement()
   {
   }


   //************* GETTING METHODS
   //************************************************************************************
   //************************************************************************************


   //************* GETTING METHODS ******************************************************
   //************************************************************************************
   //************************************************************************************


   void UpdatedLagrangianUJwPElement::GetDofList( DofsVectorType& rElementalDofList, const  ProcessInfo& rCurrentProcessInfo ) const
   {
      rElementalDofList.resize( 0 );

      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
      {
         rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
         rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );

         if( dimension == 3 )
            rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );

         rElementalDofList.push_back( GetGeometry()[i].pGetDof( JACOBIAN ) );
         rElementalDofList.push_back( GetGeometry()[i].pGetDof( WATER_PRESSURE ) );
      }
   }


   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUJwPElement::EquationIdVector( EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo ) const
   {
      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
      unsigned int element_size          = number_of_nodes * dimension + 2*number_of_nodes;

      if ( rResult.size() != element_size )
         rResult.resize( element_size, false );

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         int index = i * ( dimension + 2);
         rResult[index]     = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
         rResult[index + 1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();

         if( dimension == 3)
         {
            rResult[index + 2] = GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId();
            rResult[index + 3] = GetGeometry()[i].GetDof( JACOBIAN ).EquationId();
            rResult[index + 4] = GetGeometry()[i].GetDof( WATER_PRESSURE ).EquationId();
         }
         else
         {
            rResult[index + 2] = GetGeometry()[i].GetDof( JACOBIAN ).EquationId();
            rResult[index + 3] = GetGeometry()[i].GetDof( WATER_PRESSURE ).EquationId();
         }

      }

   }

   //*********************************DISPLACEMENT***************************************
   //************************************************************************************

   void UpdatedLagrangianUJwPElement::GetValuesVector( Vector& rValues, int Step ) const
   {
      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
      unsigned int       element_size    = number_of_nodes * dimension + 2*number_of_nodes;

      if ( rValues.size() != element_size ) rValues.resize( element_size, false );


      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int index = i * (dimension +2) ;
         rValues[index]     = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X, Step );
         rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y, Step );

         if ( dimension == 3 )
         {
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Z, Step );
            rValues[index + 3] = GetGeometry()[i].GetSolutionStepValue( JACOBIAN, Step );
            rValues[index + 4] = GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE, Step);
         }
         else
         {
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( JACOBIAN, Step );
            rValues[index + 3] = GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE, Step);
         }

      }
   }


   //************************************VELOCITY****************************************
   //************************************************************************************

   void UpdatedLagrangianUJwPElement::GetFirstDerivativesVector( Vector& rValues, int Step ) const
   {
      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
      unsigned int       element_size    = number_of_nodes * dimension + 2*number_of_nodes;

      if ( rValues.size() != element_size ) rValues.resize( element_size, false );

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int index = i * (dimension + 2);
         rValues[index]     = GetGeometry()[i].GetSolutionStepValue( VELOCITY_X, Step );
         rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Y, Step );
         if ( dimension == 3 )
         {
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Z, Step );
            rValues[index + 3] = 0;
            rValues[index + 4] = 0;
         }
         else
         {
            rValues[index + 2] = 0;
            rValues[index + 3] = 0;
         }
      }
   }

   //*********************************ACCELERATION***************************************
   //************************************************************************************

   void UpdatedLagrangianUJwPElement::GetSecondDerivativesVector( Vector& rValues, int Step ) const
   {
      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
      unsigned int       element_size    = number_of_nodes * dimension + 2*number_of_nodes;

      if ( rValues.size() != element_size ) rValues.resize( element_size, false );


      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int index = i * (dimension + 2);
         rValues[index]     = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_X, Step );
         rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Y, Step );

         if ( dimension == 3 )
         {
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Z, Step );
            rValues[index + 3] = 0;
            rValues[index + 4] = 0;
         }
         else
         {
            rValues[index + 2] = 0;
            rValues[index + 3] = 0;
         }
      }

   }



   //************************************************************************************
   //************************************************************************************

   int  UpdatedLagrangianUJwPElement::Check( const ProcessInfo& rCurrentProcessInfo ) const
   {
      KRATOS_TRY

      int correct = 0;

      correct = UpdatedLagrangianUJElement::Check(rCurrentProcessInfo);


      //verify compatibility with the constitutive law
      ConstitutiveLaw::Features LawFeatures;
      this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetLawFeatures(LawFeatures);

      if(LawFeatures.mOptions.Is(ConstitutiveLaw::U_P_LAW))
         KRATOS_THROW_ERROR( std::logic_error, "constitutive law is not compatible with the U-J element type ", " UpdatedLagrangianUJwPElement" )

            //verify that the variables are correctly initialized

            if ( WATER_PRESSURE.Key() == 0 )
               KRATOS_THROW_ERROR( std::invalid_argument, "PRESSURE has Key zero! (check if the application is correctly registered", "" )

                  return correct;

      KRATOS_CATCH( "" );
   }


   //************* STARTING - ENDING  METHODS
   //************************************************************************************
   //************************************************************************************



   void UpdatedLagrangianUJwPElement::InitializeElementData ( ElementDataType & rVariables, const ProcessInfo& rCurrentProcessInfo)
   {
      UpdatedLagrangianUJElement::InitializeElementData( rVariables, rCurrentProcessInfo );

      mTimeStep = rCurrentProcessInfo[DELTA_TIME];

      //stabilization factor
     double StabilizationFactor = 1.0;
     if( GetProperties().Has(STABILIZATION_FACTOR_WP) ){
       StabilizationFactor = GetProperties()[STABILIZATION_FACTOR_WP];
   }
     else if( rCurrentProcessInfo.Has(STABILIZATION_FACTOR_WP) ){
       StabilizationFactor = rCurrentProcessInfo[STABILIZATION_FACTOR_WP];
            }
     GetProperties().SetValue(STABILIZATION_FACTOR_WP, StabilizationFactor);

         }
   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUJwPElement::InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
         VectorType& rRightHandSideVector,
         Flags& rCalculationFlags)

   {

      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

      //resizing as needed the LHS
      unsigned int MatSize = number_of_nodes * dimension + 2*number_of_nodes;

      if ( rCalculationFlags.Is(LargeDisplacementElement::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
      {
         if ( rLeftHandSideMatrix.size1() != MatSize )
            rLeftHandSideMatrix.resize( MatSize, MatSize, false );

         noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
      }


      //resizing as needed the RHS
      if ( rCalculationFlags.Is(LargeDisplacementElement::COMPUTE_RHS_VECTOR) ) //calculation of the matrix is required
      {
         if ( rRightHandSideVector.size() != MatSize )
            rRightHandSideVector.resize( MatSize, false );

         rRightHandSideVector = ZeroVector( MatSize ); //resetting RHS
      }
   }


   //************* COMPUTING  METHODS
   //************************************************************************************
   //************************************************************************************


   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUJwPElement::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, ElementDataType& rVariables, double& rIntegrationWeight)
   {

      KRATOS_TRY

      // define some variables
      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
      rVariables.detF0 *= rVariables.detF;
      double DeterminantF = rVariables.detF;
      rVariables.detF = 1.0;

      //contributions of the stiffness matrix calculated on the reference configuration
      MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();

      // ComputeBaseClass LHS
      LocalSystemComponents UJLocalSystem;
      unsigned int MatSize = number_of_nodes * ( dimension+1);
      MatrixType  LocalLeftHandSideMatrix = ZeroMatrix(MatSize,MatSize) ;
      UJLocalSystem.SetLeftHandSideMatrix( LocalLeftHandSideMatrix);

      // LHS. base class
      UpdatedLagrangianUJElement::CalculateAndAddLHS( UJLocalSystem, rVariables, rIntegrationWeight);


      // Reshape the BaseClass LHS and Add the Hydro Part
      WaterPressureJacobianUtilities WaterUtility;
      Matrix TotalF = prod ( rVariables.F, rVariables.F0);
      int number_of_variables = dimension + 2; // displ - Jacobian - waterPressure
      Vector VolumeForce;
      VolumeForce= this->CalculateVolumeForce( VolumeForce, rVariables);
      // 1. Create (make pointers) variables
      WaterPressureUtilities::HydroMechanicalVariables HMVariables(GetGeometry(), GetProperties() );

      HMVariables.SetBMatrix( rVariables.B);
      HMVariables.SetShapeFunctionsDerivatives( rVariables.DN_DX);
      HMVariables.SetDeformationGradient( TotalF);
      HMVariables.SetVolumeForce( VolumeForce);
      HMVariables.SetShapeFunctions( rVariables.N);

      HMVariables.DeltaTime = mTimeStep;
      HMVariables.detF0 = rVariables.detF0;
      //HMVariables.CurrentRadius
      //HMVariables.ConstrainedModulus
      HMVariables.number_of_variables = number_of_variables;

      rLeftHandSideMatrix = WaterUtility.CalculateAndAddHydromechanicalLHS( HMVariables, rLeftHandSideMatrix, LocalLeftHandSideMatrix, rIntegrationWeight);

     // LHS. water pressure stabilization
      ProcessInfo SomeProcessInfo;
      std::vector<double> Mmodulus;
      this->CalculateOnIntegrationPoints( M_MODULUS, Mmodulus, SomeProcessInfo);
      double ConstrainedModulus = Mmodulus[0];
      if ( ConstrainedModulus < 1e-5)
      {
         const double& YoungModulus          = GetProperties()[YOUNG_MODULUS];
         const double& nu    = GetProperties()[POISSON_RATIO];
         ConstrainedModulus =  YoungModulus * ( 1.0-nu)/(1.0+nu) / (1.0-2.0*nu);
      }
      HMVariables.ConstrainedModulus = ConstrainedModulus;

      rLeftHandSideMatrix = WaterUtility.CalculateAndAddStabilizationLHS( HMVariables, rLeftHandSideMatrix, rIntegrationWeight);

      rVariables.detF = DeterminantF;
      rVariables.detF0 /= rVariables.detF;

      KRATOS_CATCH("")
   }

   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUJwPElement::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, ElementDataType& rVariables, Vector& rVolumeForce, double& rIntegrationWeight)
   {

      KRATOS_TRY

      // define some variables
      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
      rVariables.detF0 *= rVariables.detF;
      double DeterminantF = rVariables.detF;
      rVariables.detF = 1.0;


      //contribution of the internal and external forces
      VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector();

      // Compute Base Class RHS
      LocalSystemComponents BaseClassLocalSystem;
      Vector BaseClassRightHandSideVector = ZeroVector ( (dimension+1) * number_of_nodes );
      BaseClassLocalSystem.SetRightHandSideVector(BaseClassRightHandSideVector );
      Vector VolumeForce = rVolumeForce;
      VolumeForce *= 0.0;
      UpdatedLagrangianUJElement::CalculateAndAddRHS( BaseClassLocalSystem, rVariables, VolumeForce, rIntegrationWeight);

      // Reshape the BaseClass RHS and Add the Hydro Part
      WaterPressureJacobianUtilities WaterUtility;
      Matrix TotalF = prod( rVariables.F, rVariables.F0);
      int number_of_variables = dimension+2; // displacement - Jacobian - waterPressure

      // 1. Create (make pointers) variables
      WaterPressureUtilities::HydroMechanicalVariables HMVariables(GetGeometry(), GetProperties() );

      HMVariables.SetBMatrix( rVariables.B);
      HMVariables.SetShapeFunctionsDerivatives( rVariables.DN_DX);
      HMVariables.SetDeformationGradient( TotalF);
      HMVariables.SetVolumeForce( rVolumeForce);
      HMVariables.SetShapeFunctions( rVariables.N);

      HMVariables.DeltaTime = mTimeStep;
      HMVariables.detF0 = rVariables.detF0;
      //HMVariables.CurrentRadius
      //HMVariables.ConstrainedModulus
      HMVariables.number_of_variables = number_of_variables;

      rRightHandSideVector = WaterUtility.CalculateAndAddHydromechanicalRHS( HMVariables, rRightHandSideVector, BaseClassRightHandSideVector, rIntegrationWeight);

      // Add Stab term
      ProcessInfo SomeProcessInfo;
      std::vector< double> Mmodulus;
      this->CalculateOnIntegrationPoints( M_MODULUS, Mmodulus, SomeProcessInfo);
      double ConstrainedModulus = Mmodulus[0];
      if ( ConstrainedModulus < 1e-5)
      {
         const double& YoungModulus          = GetProperties()[YOUNG_MODULUS];
         const double& nu    = GetProperties()[POISSON_RATIO];
         ConstrainedModulus =  YoungModulus * ( 1.0-nu)/(1.0+nu) / (1.0-2.0*nu);
      }
      HMVariables.ConstrainedModulus = ConstrainedModulus;
      rRightHandSideVector = WaterUtility.CalculateAndAddStabilization( HMVariables, rRightHandSideVector, rIntegrationWeight);


      rVariables.detF = DeterminantF;
      rVariables.detF0 /= rVariables.detF;

      KRATOS_CATCH( "" )
   }


   //************************************************************************************
   //************************************************************************************
   // ** mass and damping matrix, copied directly from LargeDisplacementUPElement
   void UpdatedLagrangianUJwPElement::CalculateMassMatrix( MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo )
   {
      KRATOS_TRY

      //lumped
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
      const unsigned int number_of_nodes = GetGeometry().size();
      unsigned int MatSize = number_of_nodes * (dimension + 2); // displacement-J_wP

      if ( rMassMatrix.size1() != MatSize )
         rMassMatrix.resize( MatSize, MatSize, false );

      rMassMatrix = ZeroMatrix( MatSize, MatSize );

      // Not Lumped Mass Matrix (numerical integration):

      //reading integration points
      IntegrationMethod CurrentIntegrationMethod = mThisIntegrationMethod; //GeometryData::IntegrationMethod::GI_GAUSS_2; //GeometryData::IntegrationMethod::GI_GAUSS_1;

      const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( CurrentIntegrationMethod  );

      ElementDataType Variables;
      this->InitializeElementData(Variables,rCurrentProcessInfo);


      for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
   {
         //compute element kinematics
         this->CalculateKinematics( Variables, PointNumber );

         //getting informations for integration
         double IntegrationWeight = integration_points[PointNumber].Weight() * Variables.detJ;

         IntegrationWeight = this->CalculateIntegrationWeight( IntegrationWeight );  // multiplies by thickness

         //compute the current density
         double Determinant = Variables.detF * Variables.detF0;

      double WaterDensity = GetProperties()[DENSITY_WATER];
         double MixtureDensity = GetProperties()[DENSITY];

         double CurrentDensity = MixtureDensity / Determinant + ( Determinant-1)/Determinant*WaterDensity;

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
            unsigned int indexupi = ( dimension + 2 ) * i;

         for ( unsigned int j = 0; j < number_of_nodes; j++ )
         {
               unsigned int indexupj = ( dimension + 2 ) * j;

            for ( unsigned int k = 0; k < dimension; k++ )
            {
                  rMassMatrix( indexupi+k , indexupj+k ) += Variables.N[i] * Variables.N[j] * CurrentDensity * IntegrationWeight;
            }
            }
      }

   }


      KRATOS_CATCH("")
   }

   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUJwPElement::CalculateDampingMatrix( MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo )
   {
      KRATOS_TRY

      //0.-Initialize the DampingMatrix:
      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

      //resizing as needed the LHS
      unsigned int MatSize = number_of_nodes * ( dimension + 2);

      if ( rDampingMatrix.size1() != MatSize )
         rDampingMatrix.resize( MatSize, MatSize, false );

      noalias( rDampingMatrix ) = ZeroMatrix( MatSize, MatSize );

      //1.-Calculate StiffnessMatrix:

      MatrixType LHSMatrix  = Matrix();

      this->CalculateLeftHandSide( LHSMatrix, rCurrentProcessInfo );

      MatrixType StiffnessMatrix  = Matrix();

      if ( StiffnessMatrix.size1() != MatSize )
         StiffnessMatrix.resize( MatSize, MatSize, false );

      StiffnessMatrix = ZeroMatrix( MatSize, MatSize );

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         //unsigned int indexup = i * ( dimension + 2);

         for ( unsigned int j = 0; j < dimension; j++ )
         {
            //StiffnessMatrix( indexup+j , indexup+j ) = LHSMatrix( indexup+j , indexup+j ); // això és una mica rotllo perquè agafa Kuu, que és diferent de K de displacement-based.
         }
      }

      //2.-Calculate MassMatrix:

      MatrixType MassMatrix  = Matrix();

      this->CalculateMassMatrix ( MassMatrix, rCurrentProcessInfo );


      //3.-Get Damping Coeffitients (RAYLEIGH_ALPHA, RAYLEIGH_BETA)
      double alpha = 0;
      if( GetProperties().Has(RAYLEIGH_ALPHA) ){
         alpha = GetProperties()[RAYLEIGH_ALPHA];
      }
      else if( rCurrentProcessInfo.Has(RAYLEIGH_ALPHA) ){
         alpha = rCurrentProcessInfo[RAYLEIGH_ALPHA];
      }

      double beta  = 0;
      if( GetProperties().Has(RAYLEIGH_BETA) ){
         beta = GetProperties()[RAYLEIGH_BETA];
      }
      else if( rCurrentProcessInfo.Has(RAYLEIGH_BETA) ){
         beta = rCurrentProcessInfo[RAYLEIGH_BETA];
      }

      //4.-Compose the Damping Matrix:

      //Rayleigh Damping Matrix: alpha*M + beta*K
      rDampingMatrix  = alpha * MassMatrix;
      rDampingMatrix += beta  * StiffnessMatrix;


      KRATOS_CATCH( "" )
   }

   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUJwPElement::save( Serializer& rSerializer ) const
   {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, UpdatedLagrangianUJElement )
   }

   void UpdatedLagrangianUJwPElement::load( Serializer& rSerializer )
   {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, UpdatedLagrangianUJElement )
      rSerializer.load("DeformationGradientF0",mDeformationGradientF0);
      rSerializer.load("DeterminantF0",mDeterminantF0);
   }

}
