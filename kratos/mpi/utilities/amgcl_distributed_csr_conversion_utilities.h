
//    |  /           | 
//    ' /   __| _` | __|  _ \   __| 
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/ 
//                   Multi-Physics  
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Denis Demidov
//                   Riccardo Rossi
//                   
//

#if !defined(KRATOS_DISTRIBUTED_CSR_CONVERSION_UTILITIES_H_INCLUDED)
#define  KRATOS_DISTRIBUTED_CSR_CONVERSION_UTILITIES_H_INCLUDED

#include "containers/distributed_csr_matrix.h"
#include "amgcl/backend/builtin.hpp"
#include "amgcl/backend/interface.hpp"
#include "amgcl/mpi/distributed_matrix.hpp"
#include "amgcl/adapter/zero_copy.hpp"
#include "utilities/amgcl_csr_conversion_utilities.h"

namespace Kratos
{

/**
Utilities to convert the distributed_csr matrix to other libraries
 */
class AmgclDistributedCSRConversionUtilities
{

public:

    /**
     This function returns a shared pointer to an Amgcl distributed_matrix
     */
	template< class TDataType, class TIndexType >
	static Kratos::shared_ptr<amgcl::mpi::distributed_matrix<amgcl::backend::builtin<double>>> ConvertToAmgcl(
			const DistributedCsrMatrix<TDataType, TIndexType>& rA,
			DenseVector<TIndexType>& global_index2,
			bool MoveToBackend=true)
	{
		IndexType chunk = rA.local_size1();

    	auto loc_a = amgcl::adapter::zero_copy(chunk, 
                rA.GetDiagonalBlock().index1_data().begin(),
                rA.GetDiagonalBlock().index2_data().begin(),
                rA.GetDiagonalBlock().value_data().begin()
			);
		loc_a->ncols = rA.GetDiagonalBlock().size2(); //important if the matrix is not square

		auto rem_a = amgcl::adapter::zero_copy(chunk, 
			rA.GetOffDiagonalBlock().index1_data().begin(),
			global_index2.data().begin(),
			rA.GetOffDiagonalBlock().value_data().begin()
		);
		loc_a->ncols = rA.GetDiagonalBlock().size2(); //important if the matrix is not square

		auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator( rA.GetComm());
		amgcl::mpi::communicator comm(raw_mpi_comm);

		auto pAmgcl = Kratos::make_shared<amgcl::mpi::distributed_matrix<amgcl::backend::builtin<double>>>(comm, loc_a, rem_a);

		if(MoveToBackend)
			pAmgcl->move_to_backend();

		return pAmgcl;
	}

	template< class TDataType, class TIndexType >
	static DistributedCsrMatrix<TDataType, TIndexType> ConvertToCsrMatrix(
			amgcl::mpi::distributed_matrix<amgcl::backend::builtin<double>>& rA //cannot be made const since i need to modify some data in-place
			)
	{
		if(!rA.local())
			KRATOS_ERROR << "matrix A was moved to backend, so it is impossible to convert it back to CSR matrix" << std::endl;

		DistributedCsrMatrix<TDataType, TIndexType> Aconverted;

		DataCommunicator& kratos_comm=ParallelEnvironment::GetDefaultDataCommunicator();

		rA.local()->own_data=false;
		rA.remote()->own_data=false;

		//create row numbering and col numbering
		Aconverted.pGetRowNumbering() = Kratos::make_unique<DistributedNumbering<IndexType>>(kratos_comm,rA.local()->nrows);
		Aconverted.pGetColNumbering() = Kratos::make_unique<DistributedNumbering<IndexType>>(kratos_comm,rA.local()->ncols);

		//here we fill the global to local mapping
		const auto& comm_pattern = rA.cpat();
		for(TIndexType i=0; i<rA.remote()->nnz; ++i){  //TODO: i suspect there is a smarter way to do this
			TIndexType id = rA.remote()->col[i]; 
			TIndexType local_id = comm_pattern.local_index(id);
			Aconverted.GetOffDiagonalLocalIds()[id] = local_id; 
			rA.remote()->col[i] = local_id; //note that here we overwrite the amgcl data!
		}

		//here we fill the local to global mapping (the inverse of the previous one)
		Aconverted.GetOffDiagonalGlobalIds().resize(Aconverted.GetOffDiagonalLocalIds().size());

		for(auto item : Aconverted.GetOffDiagonalLocalIds()){
			Aconverted.GetOffDiagonalGlobalIds()[item.second] = item.first; //item.second=local_id, item.first=global_id
		}

		//setting col size for both matrices
		Aconverted.GetDiagonalBlock().SetColSize(rA.local()->ncols); 
		Aconverted.GetOffDiagonalBlock().SetColSize(Aconverted.GetOffDiagonalGlobalIds().size()); 

		//convert diagonal block
		Aconverted.GetDiagonalBlock() = std::move(AmgclCSRConversionUtilities::ConvertToCsrMatrix<TDataType,TIndexType>(*(rA.local().get())));

		//convert off diagonal block. Note that we need to change the usage of the index2. We do it in place
		Aconverted.GetOffDiagonalBlock() = std::move(AmgclCSRConversionUtilities::ConvertToCsrMatrix<TDataType,TIndexType>(*(rA.remote().get())));

		//fill the vector importer, so that the matrix can be used to do calculations
		Aconverted.pGetVectorImporter() = Kratos::make_unique<DistributedVectorImporter<TDataType,IndexType>>(
			kratos_comm,Aconverted.GetOffDiagonalGlobalIds(), 
			Aconverted.GetColNumbering()
			); 

		//the following data are simply not available in Amgcl.
		//warning: THE IMPORTANT IMPLICATION IS THAT THE MATRIX CANNOT BE USED FOR ASSEMBLY WITHOUT PROVIDING MORE INFORMATION
		// mSendCachedIJ
        // mRecvCachedIJ
        // mOffDiagonalLocalIds
        // mOffDiagonalGlobalIds
        // mfem_assemble_colors

		//set ownership
        if(rA.local()->own_data == false){ //if rA is not the owner, Aconverted cannot be
            Aconverted.GetDiagonalBlock().SetIsOwnerOfData(false);
        }
        else{ //if rA is the owner, transfer ownership to the csr_matrix
            rA.local()->own_data = false;
            Aconverted.GetDiagonalBlock().SetIsOwnerOfData(true);
        }		

        if(rA.remote()->own_data == false){ //if rA is not the owner, Aconverted cannot be
            Aconverted.GetOffDiagonalBlock().SetIsOwnerOfData(false);
        }
        else{ //if rA is the owner, transfer ownership to the csr_matrix
            rA.local()->own_data = false;
            Aconverted.GetOffDiagonalBlock().SetIsOwnerOfData(true);
        }		

		return Aconverted;


	}	

};

}

#endif // KRATOS_DISTRIBUTED_CSR_CONVERSION_UTILITIES_H_INCLUDED