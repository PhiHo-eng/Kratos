/*
==============================================================================
KratosULFApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Pawel Ryzhakov
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/


//
//   Project Name:        Kratos
//   Last Modified by:    $Author: anonymous $
//   Date:                $Date: 2007-03-28 10:38:34 $
//   Revision:            $Revision: 1.1 $
//
//


#if !defined(KRATOS_UPDATED_LAGRANGIAN__FLUID3D_H_INCLUDED )
#define  KRATOS_UPDATED_LAGRANGIAN__FLUID3D_H_INCLUDED



// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/serializer.h"

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

/// Short class definition.
/** Detail class definition.
*/
class UpdatedLagrangianFluid3D
    : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of UpdatedLagrangianFluid3D
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(UpdatedLagrangianFluid3D);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    UpdatedLagrangianFluid3D(IndexType NewId, GeometryType::Pointer pGeometry);
    UpdatedLagrangianFluid3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

    /// Destructor.
    ~UpdatedLagrangianFluid3D() override;


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
    //virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);

    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

    void GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo);

    void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo);

    void Calculate(const Variable<double >& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo);
    void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo);
    void CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo);
    void GetValuesVector(Vector& values, int Step = 0);
    void GetFirstDerivativesVector(Vector& values, int Step = 0);
    void GetSecondDerivativesVector(Vector& values, int Step = 0);


    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "UpdatedLagrangianFluid3D #" ;
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info() << Id();
    }

    /// Print object's data.
//      virtual void PrintData(std::ostream& rOStream) const;


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
    /*
    static boost::numeric::ublas::bounded_matrix<double,4,4> msMassFactors;
    static boost::numeric::ublas::bounded_matrix<double,4,3> msDN_Dx;
    static array_1d<double,4> msN; //dimension = number of nodes
    //static Matrix msDN_DX;
    //static Matrix msMassFactors;
    static array_1d<double,3> ms_vel_gauss; //dimesion coincides with space dimension
    static array_1d<double,4> ms_temp_vec_np; //dimension = number of nodes
    static array_1d<double,4> ms_u_DN;

    static boost::numeric::ublas::bounded_matrix<double,6,12> msB;
    static boost::numeric::ublas::bounded_matrix<double,6,6> ms_constitutive_matrix;
    static boost::numeric::ublas::bounded_matrix<double,6,12> ms_temp;

    static array_1d<double,6> ms_temp_vec;
    */


    ///@}
    ///@name Member Variables
    ///@{
    double mA0; //stores area (or Volume in 3D) at the beginning of time-step

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    // A private default constructor necessary for serialization
    UpdatedLagrangianFluid3D() : Element()
    {
    }

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    }



    ///@}
    ///@name Private Operators
    ///@{
    void Stage1(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, unsigned int ComponentIndex);
    void Stage2(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

    void MeshMovingStep(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, unsigned int ComponentIndex);

    //inline void CalculateGeometryData(Matrix& msDN_DX, Vector& N, double& Area)
    inline void CalculateGeometryData(BoundedMatrix<double,3,2>& DN_DX, array_1d<double,3>& N, double& Area);

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

    /// Assignment operator.
    //UpdatedLagrangianFluid3D& operator=(const UpdatedLagrangianFluid3D& rOther);

    /// Copy constructor.
    //UpdatedLagrangianFluid3D(const UpdatedLagrangianFluid3D& rOther);


    ///@}

}; // Class UpdatedLagrangianFluid3D

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
				    UpdatedLagrangianFluid3D& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
				    const UpdatedLagrangianFluid3D& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

}  // namespace Kratos.

#endif // KRATOS_UPDATED_LAGRANGIAN__FLUID3D_H_INCLUDED   defined 


