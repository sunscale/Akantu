/**
 * @file   petsc_matrix.cc
 * @author Aurelia Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Mon Jul 21 17:40:41 2014
 *
 * @brief  Implementation of PETSc matrix class
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "petsc_matrix.hh"
#include "static_communicator.hh"
#include "static_communicator_mpi.hh"
#include "mpi_type_wrapper.hh"
#include "dof_synchronizer.hh"
/* -------------------------------------------------------------------------- */
#include <cstring>
#include <mpi.h>
#include <petscmat.h>
#include <petscao.h>
#include <petscis.h>
#include <petscsys.h>  


__BEGIN_AKANTU__

struct PETScWrapper {
  Mat mat;
  AO ao;
  ISLocalToGlobalMapping mapping;
  /// pointer to the MPI communicator for PETSc commands
  MPI_Comm communicator;
};

/* -------------------------------------------------------------------------- */
PETScMatrix::PETScMatrix(UInt size,
			 const SparseMatrixType & sparse_matrix_type,
			 const ID & id,
			 const MemoryID & memory_id) :
  SparseMatrix(size, sparse_matrix_type, id, memory_id), 
  petsc_wrapper(new PETScWrapper), d_nnz(0,1,"dnnz"), o_nnz(0,1,"onnz") {
  AKANTU_DEBUG_IN();
  
  this->init();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
PETScMatrix::PETScMatrix(const SparseMatrix & matrix,
			   const ID & id,
			   const MemoryID & memory_id) :
  SparseMatrix(matrix, id, memory_id),
  petsc_wrapper(new PETScWrapper), d_nnz(0,1,"dnnz"), o_nnz(0,1,"onnz") {
  AKANTU_DEBUG_IN();

  this->init();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
PETScMatrix::~PETScMatrix() {
  AKANTU_DEBUG_IN();

  PetscErrorCode ierr;

  ierr = MatDestroy(&(this->petsc_wrapper->mat)); CHKERRXX(ierr);
  delete petsc_wrapper;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void PETScMatrix::init() {
  AKANTU_DEBUG_IN();

#if defined(AKANTU_USE_MPI)
    StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
    const StaticCommunicatorMPI & mpi_st_comm =
      dynamic_cast<const StaticCommunicatorMPI &>(comm.getRealStaticCommunicator());
    this->petsc_wrapper->communicator = mpi_st_comm.getMPITypeWrapper().getMPICommunicator();
#else
    this->petsc_wrapper->communicator = PETSC_COMM_SELF;
#endif
 
  PetscErrorCode ierr;

  ierr = MatCreate(this->petsc_wrapper->communicator, &(this->petsc_wrapper->mat)); CHKERRXX(ierr);
  /// set matrix type
  ierr = MatSetType(this->petsc_wrapper->mat,  MATAIJ); CHKERRXX(ierr);


  if (this->sparse_matrix_type==_symmetric)
    /// set flag for symmetry to enable ICC/Cholesky preconditioner
    ierr = MatSetOption(this->petsc_wrapper->mat, MAT_SYMMETRIC, PETSC_TRUE); CHKERRXX(ierr);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void PETScMatrix::resize(const DOFSynchronizer & dof_synchronizer) {
  AKANTU_DEBUG_IN();

  PetscErrorCode ierr;

  this->dof_synchronizer = &const_cast<DOFSynchronizer &>(dof_synchronizer);

  /// find the number of dofs corresponding to master or local nodes
  UInt nb_dofs = dof_synchronizer.getNbDOFs();
  UInt nb_local_master_dofs = 0;

  /// create array to store the global equation number of all local and master dofs
  Array<Int> local_master_eq_nbs(nb_dofs);
  Array<Int>::scalar_iterator it_eq_nb = local_master_eq_nbs.begin();

  Array<Int> dof_types = dof_synchronizer.getDOFTypes();
  Array<Int>::scalar_iterator it_dof_type = dof_types.begin();

  /// get the pointer to the gloabl equation number array
  Int * eq_nb_val = dof_synchronizer.getGlobalDOFEquationNumbers().storage();

  for (UInt i = 0; i <nb_dofs; ++i, ++it_dof_type) {
    if (*it_dof_type == -2 || *it_dof_type == -1 ) {
      /// get for the local dof the global equation number
      *it_eq_nb = eq_nb_val[i];
      ++it_eq_nb;
      ++nb_local_master_dofs;
    }
  }

  local_master_eq_nbs.resize(nb_local_master_dofs);

  /// set the local size
  this->local_size = nb_local_master_dofs;

  /// resize PETSc matrix
#if defined(AKANTU_USE_MPI)
  ierr = MatSetSizes(this->petsc_wrapper->mat, this->size, this->size, this->local_size, this->local_size); CHKERRXX(ierr);
#else
  ierr =  MatSetSizes(this->petsc_wrapper->mat, this->local_size, this->local_size); CHKERRXX(ierr);
#endif
  

  /// create mapping from akantu global numbering to petsc global numbering
  this->createGlobalAkantuToPETScMap(local_master_eq_nbs.storage());
 
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void PETScMatrix::createGlobalAkantuToPETScMap(Int* local_master_eq_nbs_ptr) {
  AKANTU_DEBUG_IN();

  PetscErrorCode ierr;
 
  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator(); 
  UInt rank = comm.whoAmI();

  //initialize vector to store the number of local and master nodes that are assigned to each processor
  Vector<UInt> master_local_ndofs_per_proc(nb_proc);

  /// store the nb of master and local dofs on each processor
  master_local_ndofs_per_proc(rank) = this->local_size;
  


  /// exchange the information among all processors
  comm.allGather(master_local_ndofs_per_proc.storage(), 1);

  /// each processor creates a map for his akantu global dofs to the corresponding petsc global dofs

  /// determine the PETSc-index for the first dof on each processor
  UInt petsc_start_index = 0;

  for (UInt i = 0; i < rank; ++i) {
    petsc_start_index +=  master_local_ndofs_per_proc(rank);
  }

  /// create array for petsc ordering
  Array<Int> petsc_dofs(this->local_size);
  Array<Int>::scalar_iterator it_petsc_dofs = petsc_dofs.begin();

  for (UInt i = petsc_start_index; i < petsc_start_index + this->local_size; ++i, ++it_petsc_dofs) {
    *it_petsc_dofs = i; 
  }

  ierr = AOCreateBasic(this->petsc_wrapper->communicator, this->local_size, local_master_eq_nbs_ptr, petsc_dofs.storage(), &(this->petsc_wrapper->ao)); CHKERRXX(ierr);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
// void PETScMatrix::createLocalAkantuToPETScMap(const DOFSynchronizer & dof_synchronizer) {
//   AKANTU_DEBUG_IN();

//   AKANTU_DEBUG_ASSERT(this->petsc_wrapper->ao != NULL,
//                       "You should first create a mapping from the global"
// 		      << " Akantu numbering to the global PETSc numbering");

//   this->dof_synchronizer = &const_cast<DOFSynchronizer &>(dof_synchronizer);

//   /// get the number of dofs
//   Int nb_dofs = dof_synchronizer.getNbDOFs();

//   /// get the global equation numbers for each node
//   Array<Int> global_dof_equation_numbers = dof_synchronizer.getGlobalDOFEquationNumbers();

//   /// map the global dof equation numbers to the corresponding PETSc ordering
//   AOApplicationToPETSc(this->petsc_wrapper->ao, nb_dofs,
// 		       global_dof_equation_numbers.storage());

//   /// create the mapping from the local Akantu ordering to the global PETSc ordering
//   ISLocalToGlobalMappingCreate(this->petsc_wrapper->communicator,
// 			       1, nb_dofs, global_dof_equation_numbers.storage(),
// 			       PETSC_COPY_VALUES, mapping);

//   AKANTU_DEBUG_OUT();
// }

/* -------------------------------------------------------------------------- */
void PETScMatrix::buildProfile(const Mesh & mesh, const DOFSynchronizer & dof_synchronizer, UInt nb_degree_of_freedom) {
  AKANTU_DEBUG_IN();

  //clearProfile();

  PetscErrorCode ierr;

  /// resize arrays to store the number of nonzeros in each row
  this->d_nnz.resize(this->local_size);
  this->o_nnz.resize(this->local_size);

  /// set arrays to zero everywhere
  this->d_nnz.set(0);
  this->o_nnz.set(0); 

  this->dof_synchronizer = &const_cast<DOFSynchronizer &>(dof_synchronizer);


  /// get the types for each local dof
  const Array<Int> & dof_types = dof_synchronizer.getDOFTypes();
  Int* dof_type_ptr = dof_types.storage();
  // if(irn_jcn_to_k) delete irn_jcn_to_k;
  // irn_jcn_to_k = new std::map<std::pair<UInt, UInt>, UInt>;
  

  this->dof_synchronizer = &const_cast<DOFSynchronizer &>(dof_synchronizer);

  coordinate_list_map::iterator irn_jcn_k_it;

  Int * eq_nb_val = dof_synchronizer.getGlobalDOFEquationNumbers().storage();

  /// Loop over all the ghost types
  for (ghost_type_t::iterator gt = ghost_type_t::begin(); gt != ghost_type_t::end(); ++gt) {
    const GhostType & ghost_type = *gt;
    /// Loop over all element types 
    Mesh::type_iterator it  = mesh.firstType(mesh.getSpatialDimension(), ghost_type, _ek_not_defined);
    Mesh::type_iterator end = mesh.lastType (mesh.getSpatialDimension(), ghost_type, _ek_not_defined);

    for(; it != end; ++it) {
      UInt nb_element = mesh.getNbElement(*it);
      UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(*it);
      UInt size_mat = nb_nodes_per_element * nb_degree_of_freedom;

      UInt * conn_val = mesh.getConnectivity(*it, _not_ghost).storage();
      Int * global_eq_nb_val = new Int[nb_degree_of_freedom * nb_nodes_per_element];
      Int * local_eq_nb_val = new Int[nb_degree_of_freedom * nb_nodes_per_element];


      for (UInt e = 0; e < nb_element; ++e) {
	Int * tmp_global_eq_nb_val = global_eq_nb_val;
	Int * tmp_local_eq_nb_val = local_eq_nb_val;

	for (UInt i = 0; i < nb_nodes_per_element; ++i) {
	  UInt n = conn_val[i];
	  for (UInt d = 0; d < nb_degree_of_freedom; ++d) {
	    *tmp_global_eq_nb_val++ = eq_nb_val[n * nb_degree_of_freedom + d];
	    *tmp_local_eq_nb_val++ = n * nb_degree_of_freedom + d;
	  }
	  // memcpy(tmp_global_eq_nb_val, eq_nb_val + n * nb_degree_of_freedom, nb_degree_of_freedom * sizeof(Int));
	  // tmp_global_eq_nb_val += nb_degree_of_freedom;
	}

	for (UInt i = 0; i < size_mat; ++i) {
	  if (dof_type_ptr[local_eq_nb_val[i]] == -2 || dof_type_ptr[local_eq_nb_val[i]] == -1 ) {
	    UInt c_irn = global_eq_nb_val[i];
	    if(c_irn < size ) {
	      //UInt j_start = (sparse_matrix_type == _symmetric) ? i : 0;
	      for (UInt j = 0; j < size_mat; ++j) {
		UInt c_jcn = global_eq_nb_val[j];
		if(c_jcn < size) {
		  KeyCOO irn_jcn = key(c_irn, c_jcn);
		  irn_jcn_k_it = irn_jcn_k.find(irn_jcn);

		  if (irn_jcn_k_it == irn_jcn_k.end()) {
		    if (dof_type_ptr[local_eq_nb_val[j]] == -2 || dof_type_ptr[local_eq_nb_val[j]] == -1 )
		      this->d_nnz(local_eq_nb_val[i]) += 1;
		    else
		      this->o_nnz(local_eq_nb_val[i]) += 1;
		  }
		}
	      }
	    }
	  }
	}
	conn_val += nb_nodes_per_element;
      }

      delete [] global_eq_nb_val;
      delete [] local_eq_nb_val;
    }

    // /// for pbc @todo correct it for parallel
    // if(StaticCommunicator::getStaticCommunicator().getNbProc() == 1) {
    //   for (UInt i = 0; i < size; ++i) {
    // 	KeyCOO irn_jcn = key(i, i);
    // 	irn_jcn_k_it = irn_jcn_k.find(irn_jcn);
    // 	if(irn_jcn_k_it == irn_jcn_k.end()) {
    // 	  irn_jcn_k[irn_jcn] = nb_non_zero;
    // 	  irn.push_back(i + 1);
    // 	  jcn.push_back(i + 1);
    // 	  nb_non_zero++;
    // 	}
    //   }
    // }
  }

  /// map arrays of nonzeros per row from akantu global to petsc global numbering
  ierr = AOApplicationToPetsc(this->petsc_wrapper->ao,
		       local_size, d_nnz.storage()); CHKERRXX(ierr);
  ierr = AOApplicationToPetsc(this->petsc_wrapper->ao,
		       local_size, o_nnz.storage()); CHKERRXX(ierr);

  // std::string mat_type;
  // mat_type.resize(20, 'a');
  // std::cout << "MatType: " << mat_type << std::endl;
  // const char * mat_type_ptr = mat_type.c_str();

  MatType type;
  MatGetType(this->petsc_wrapper->mat, &type);

  //  PetscTypeCompare((PetscObject)(this->petsc_wrapper->mat),MATSEQAIJ,&sametype);
  //ierr = PetscTypeCompare(pet, MATSEQAIJ, &sametype); CHKERRXX(ierr);

  /// build profile:
  if (strcmp(type, MATSEQAIJ) == 0) {
    ierr = MatSeqAIJSetPreallocation(this->petsc_wrapper->mat,
				     0, d_nnz.storage()); CHKERRXX(ierr);
  } else if ((strcmp(type, MATMPIAIJ) == 0)) {
    ierr = MatMPIAIJSetPreallocation(this->petsc_wrapper->mat,
				     0, d_nnz.storage(), 0,
				     o_nnz.storage()); CHKERRXX(ierr); 
  } else {
    AKANTU_DEBUG_ERROR("The type " << type << " of PETSc matrix is not handled by"
		       << " akantu in the preallocation step");
  }

  /// assemble matrix
  MatAssemblyBegin(this->petsc_wrapper->mat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(this->petsc_wrapper->mat, MAT_FINAL_ASSEMBLY);

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void PETScMatrix::saveMatrix(const std::string & filename) const {
  AKANTU_DEBUG_IN();

  PetscErrorCode ierr;
  
  /// create Petsc viewer
  PetscViewer viewer; 
  ierr = PetscViewerASCIIOpen(this->petsc_wrapper->communicator, filename.c_str(), &viewer); CHKERRXX(ierr);
  /// save the matrix
  ierr =  MatView(this->petsc_wrapper->mat, viewer); CHKERRXX(ierr);

  /// destroy the viewer
  ierr =  PetscViewerDestroy(&viewer); CHKERRXX(ierr);

  AKANTU_DEBUG_OUT();
}

// /* -------------------------------------------------------------------------- */
// Array<Real> & operator*=(Array<Real> & vect, const PETScMatrix & mat) {
//   AKANTU_DEBUG_IN();

//   Array<Real> result 
//   MatMult(this->mat, vect, )

//   // AKANTU_DEBUG_ASSERT((vect.getSize()*vect.getNbComponent() == mat.getSize()) &&
//   // 		      (vect.getNbComponent() == mat.getNbDegreOfFreedom()),
//   // 		      "The size of the matrix and the vector do not match");

//   const PETScMatrixType & sparse_matrix_type = mat.getPETScMatrixType();
//   DOFSynchronizer * dof_synchronizer = mat.getDOFSynchronizerPointer();

//   UInt nb_non_zero = mat.getNbNonZero();
//   Real * tmp = new Real [vect.getNbComponent() * vect.getSize()];
//   std::fill_n(tmp, vect.getNbComponent() * vect.getSize(), 0);

//   Int * i_val  = mat.getIRN().storage();
//   Int * j_val  = mat.getJCN().storage();
//   Real * a_val = mat.getA().storage();

//   Real * vect_val = vect.storage();

//   for (UInt k = 0; k < nb_non_zero; ++k) {
//     UInt i = *(i_val++);
//     UInt j = *(j_val++);
//     Real a = *(a_val++);

//     UInt local_i = i - 1;
//     UInt local_j = j - 1;
//     if(dof_synchronizer) {
//       local_i = dof_synchronizer->getDOFLocalID(local_i);
//       local_j = dof_synchronizer->getDOFLocalID(local_j);
//     }

//     tmp[local_i] += a * vect_val[local_j];
//     if((sparse_matrix_type == _symmetric) && (local_i != local_j))
//       tmp[local_j] += a * vect_val[local_i];
//   }

//   memcpy(vect_val, tmp, vect.getNbComponent() * vect.getSize() * sizeof(Real));
//   delete [] tmp;

//   if(dof_synchronizer)
//     dof_synchronizer->reduceSynchronize<AddOperation>(vect);

//   AKANTU_DEBUG_OUT();

//   return vect;
// }

// /* -------------------------------------------------------------------------- */
// void PETScMatrix::copyContent(const PETScMatrix & matrix) {
//   AKANTU_DEBUG_IN();
//   AKANTU_DEBUG_ASSERT(nb_non_zero == matrix.getNbNonZero(),
// 		      "The two matrices don't have the same profiles");

//   MatCopy(this->mat, matrix->mat, SAME_NONZERO_PATTERN);

//   AKANTU_DEBUG_OUT();
// }

// ///* -------------------------------------------------------------------------- */
// //void PETScMatrix::copyProfile(const PETScMatrix & matrix) {
// //  AKANTU_DEBUG_IN();
// //  irn = matrix.irn;
// //  jcn = matrix.jcn;
// //  nb_non_zero = matrix.nb_non_zero;
// //  irn_jcn_k = matrix.irn_jcn_k;
// //  a.resize(nb_non_zero);
// //  AKANTU_DEBUG_OUT();
// //}

/* -------------------------------------------------------------------------- */
void PETScMatrix::add(const SparseMatrix & matrix, Real alpha) {
  PetscErrorCode ierr;
  AKANTU_DEBUG_ASSERT(nb_non_zero == matrix.getNbNonZero(),
		      "The two matrix don't have the same profiles");

  const Array<Real> matrix_a = matrix.getA();
  Array<Int> matrix_irn = matrix.getIRN();
  Array<Int> matrix_jcn = matrix.getJCN();
  UInt matrix_nb_non_zero = matrix.getNbNonZero();

  /// map global akantu to global petsc
  ierr = AOApplicationToPetsc(this->petsc_wrapper->ao, matrix_nb_non_zero,matrix_irn.storage()); CHKERRXX(ierr);
  ierr = AOApplicationToPetsc(this->petsc_wrapper->ao, matrix_nb_non_zero,matrix_jcn.storage()); CHKERRXX(ierr);

  Real val_to_add = 0;

  for (UInt n = 0; n < matrix_nb_non_zero; ++n) {
    val_to_add = alpha * matrix_a(n);
    ierr = MatSetValue(this->petsc_wrapper->mat, matrix_irn(n), matrix_jcn(n), val_to_add,INSERT_VALUES); CHKERRXX(ierr);
  }


}

/* -------------------------------------------------------------------------- */
void PETScMatrix::performAssembly() {
  PetscErrorCode ierr;
  ierr = MatAssemblyBegin(this->petsc_wrapper->mat, MAT_FINAL_ASSEMBLY); CHKERRXX(ierr);
  ierr = MatAssemblyEnd(this->petsc_wrapper->mat, MAT_FINAL_ASSEMBLY); CHKERRXX(ierr);
}


__END_AKANTU__
