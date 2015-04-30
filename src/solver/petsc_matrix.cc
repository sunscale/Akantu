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
#include "petsc_wrapper.hh"
/* -------------------------------------------------------------------------- */
#include <cstring>
#include <petscsys.h>  

__BEGIN_AKANTU__

// struct PETScWrapper {
//   Mat mat;
//   AO ao;
//   ISLocalToGlobalMapping mapping;
//   /// pointer to the MPI communicator for PETSc commands
//   MPI_Comm communicator;
// };

/* -------------------------------------------------------------------------- */
PETScMatrix::PETScMatrix(UInt size,
			 const SparseMatrixType & sparse_matrix_type,
			 const ID & id,
			 const MemoryID & memory_id) :
  SparseMatrix(size, sparse_matrix_type, id, memory_id), 
  petsc_matrix_wrapper(new PETScMatrixWrapper),
  d_nnz(0,1,"dnnz"),
  o_nnz(0,1,"onnz"),
  first_global_index(0),
  is_petsc_matrix_initialized(false) {
  AKANTU_DEBUG_IN();
  StaticSolver::getStaticSolver().registerEventHandler(*this);
  this->offset = 0;
  this->init();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
PETScMatrix::PETScMatrix(const SparseMatrix & matrix,
			 const ID & id,
			 const MemoryID & memory_id) :
  SparseMatrix(matrix, id, memory_id),
  petsc_matrix_wrapper(new PETScMatrixWrapper), 
  d_nnz(0,1,"dnnz"), 
  o_nnz(0,1,"onnz"), 
  first_global_index(0),
  is_petsc_matrix_initialized(false) {
  // AKANTU_DEBUG_IN();
  // StaticSolver::getStaticSolver().registerEventHandler(*this);
  // this->offset = 0;
  // this->init();

  // AKANTU_DEBUG_OUT();
  AKANTU_DEBUG_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
PETScMatrix::~PETScMatrix() {
  AKANTU_DEBUG_IN();

  /// destroy all the PETSc data structures used for this matrix
  this->destroyInternalData();
  StaticSolver::getStaticSolver().unregisterEventHandler(*this);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void PETScMatrix::init() {
  AKANTU_DEBUG_IN();

  /// set the communicator that should be used by PETSc
#if defined(AKANTU_USE_MPI)
    StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
    const StaticCommunicatorMPI & mpi_st_comm =
      dynamic_cast<const StaticCommunicatorMPI &>(comm.getRealStaticCommunicator());
    this->petsc_matrix_wrapper->communicator = mpi_st_comm.getMPITypeWrapper().getMPICommunicator();
#else
    this->petsc_matrix_wrapper->communicator = PETSC_COMM_SELF;
#endif
 
  PetscErrorCode ierr;
  
  /// create the PETSc matrix object
  ierr = MatCreate(this->petsc_matrix_wrapper->communicator, &(this->petsc_matrix_wrapper->mat)); CHKERRXX(ierr);


  /**
   * Set the matrix type
   * @todo PETSc does currently not support a straightforward way to
   * apply Dirichlet boundary conditions for MPISBAIJ
   * matrices. Therefore always the entire matrix is allocated. It
   * would be possible to use MATSBAIJ for sequential matrices in case
   * memory becomes critical. Also, block matrices would give a much
   * better performance. Modify this in the future!
   */

  ierr = MatSetType(this->petsc_matrix_wrapper->mat,  MATAIJ); CHKERRXX(ierr);
 
  this->is_petsc_matrix_initialized = true;


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * With this method each processor computes the dimensions of the
 * local matrix, i.e. the part of the global matrix it is storing.
 * @param dof_synchronizer dof synchronizer that maps the local
 * dofs to the global dofs and the equation numbers, i.e., the 
 * position at which the dof is assembled in the matrix
 */
void PETScMatrix::setSize() {
  AKANTU_DEBUG_IN();

  PetscErrorCode ierr;

  /// find the number of dofs corresponding to master or local nodes
  UInt nb_dofs = this->dof_synchronizer->getNbDOFs();
  UInt nb_local_master_dofs = 0;

  /// create array to store the global equation number of all local and master dofs
  Array<Int> local_master_eq_nbs(nb_dofs);
  Array<Int>::scalar_iterator it_eq_nb = local_master_eq_nbs.begin();

  /// get the pointer to the global equation number array
  Int * eq_nb_val = this->dof_synchronizer->getGlobalDOFEquationNumbers().storage();
  
  for (UInt i = 0; i <nb_dofs; ++i) {
    if (this->dof_synchronizer->isLocalOrMasterDOF(i) ) {
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
  ierr = MatSetSizes(this->petsc_matrix_wrapper->mat, this->local_size, this->local_size, this->size, this->size); CHKERRXX(ierr);
#else
  ierr =  MatSetSizes(this->petsc_matrix_wrapper->mat, this->local_size, this->local_size); CHKERRXX(ierr);
#endif

  /// create mapping from akantu global numbering to petsc global numbering
  this->createGlobalAkantuToPETScMap(local_master_eq_nbs.storage());


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * This method generates a mapping from the global Akantu equation
 * numbering to the global PETSc dof ordering 
 * @param local_master_eq_nbs_ptr Int pointer to the array of equation
 * numbers of all local or master dofs, i.e. the row indices of the
 * local matrix
 */
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
	
  for (UInt i = 0; i < rank; ++i) {
    this->first_global_index +=  master_local_ndofs_per_proc(i);
  }

  /// create array for petsc ordering
  Array<Int> petsc_dofs(this->local_size);
  Array<Int>::scalar_iterator it_petsc_dofs = petsc_dofs.begin();
	
  for (Int i = this->first_global_index; i < this->first_global_index + this->local_size; ++i, ++it_petsc_dofs) {
    *it_petsc_dofs = i; 
  }

  ierr = AOCreateBasic(this->petsc_matrix_wrapper->communicator, this->local_size, local_master_eq_nbs_ptr, petsc_dofs.storage(), &(this->petsc_matrix_wrapper->ao)); CHKERRXX(ierr);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
// void PETScMatrix::createLocalAkantuToPETScMap(const DOFSynchronizer & dof_synchronizer) {
//   AKANTU_DEBUG_IN();

//   AKANTU_DEBUG_ASSERT(this->petsc_matrix_wrapper->ao != NULL,
//                       "You should first create a mapping from the global"
// 		      << " Akantu numbering to the global PETSc numbering");

//   PetscErrorCode ierr;

//   this->dof_synchronizer = &const_cast<DOFSynchronizer &>(dof_synchronizer);

//   /// get the number of dofs
//   Int nb_dofs = dof_synchronizer.getNbDOFs();

//   /// get the global equation numbers for each node
//   Array<Int> global_dof_equation_numbers = dof_synchronizer.getGlobalDOFEquationNumbers();

//   /// map the global dof equation numbers to the corresponding PETSc ordering
//   ierr =  AOApplicationToPETSc(this->petsc_matrix_wrapper->ao, nb_dofs,
// 		       global_dof_equation_numbers.storage()); CHKERRXX(ierr);

//   /// create the mapping from the local Akantu ordering to the global PETSc ordering
//   ierr = ISLocalToGlobalMappingCreate(this->petsc_matrix_wrapper->communicator,
// 			       1, nb_dofs, global_dof_equation_numbers.storage(),
// 			       PETSC_COPY_VALUES, &(this->petsc_matrix_wrapper-mapping)); CHKERRXX(ierr);

//   AKANTU_DEBUG_OUT();
// }

/* -------------------------------------------------------------------------- */
/**
 * This function creates the non-zero pattern of the PETSc matrix. In
 * PETSc the parallel matrix is partitioned across processors such
 * that the first m0 rows belong to process 0, the next m1 rows belong
 * to process 1, the next m2 rows belong to process 2 etc.. where
 * m0,m1,m2,.. are the input parameter 'm'. i.e each processor stores
 * values corresponding to [m x N] submatrix
 * (http://www.mcs.anl.gov/petsc/).
 * @param mesh mesh discretizing the domain we want to analyze
 * @param dof_synchronizer dof synchronizer that maps the local 
 * dofs to the global dofs and the equation numbers, i.e., the 
 * position at which the dof is assembled in the matrix
 */

void PETScMatrix::buildProfile(const Mesh & mesh, const DOFSynchronizer & dof_synchronizer, UInt nb_degree_of_freedom) {
  AKANTU_DEBUG_IN();

  //clearProfile();
  this->dof_synchronizer = &const_cast<DOFSynchronizer &>(dof_synchronizer);
  this->setSize();
  PetscErrorCode ierr;

  /// resize arrays to store the number of nonzeros in each row
  this->d_nnz.resize(this->local_size);
  this->o_nnz.resize(this->local_size);

  /// set arrays to zero everywhere
  this->d_nnz.set(0);
  this->o_nnz.set(0); 


  // if(irn_jcn_to_k) delete irn_jcn_to_k;
  // irn_jcn_to_k = new std::map<std::pair<UInt, UInt>, UInt>;
  
  coordinate_list_map::iterator irn_jcn_k_it;

  Int * eq_nb_val = dof_synchronizer.getGlobalDOFEquationNumbers().storage();
  UInt nb_global_dofs = dof_synchronizer.getNbGlobalDOFs();

  /// Loop over all the ghost types
  for (ghost_type_t::iterator gt = ghost_type_t::begin(); gt != ghost_type_t::end(); ++gt) {
    const GhostType & ghost_type = *gt;
    Mesh::type_iterator it  = mesh.firstType(mesh.getSpatialDimension(), ghost_type, _ek_not_defined);
    Mesh::type_iterator end = mesh.lastType (mesh.getSpatialDimension(), ghost_type, _ek_not_defined);
    for(; it != end; ++it) {
      UInt nb_element = mesh.getNbElement(*it, ghost_type);
      UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(*it);
      UInt size_mat = nb_nodes_per_element * nb_degree_of_freedom;

      UInt * conn_val = mesh.getConnectivity(*it, ghost_type).storage();
      Int * local_eq_nb_val = new Int[nb_degree_of_freedom * nb_nodes_per_element];


      for (UInt e = 0; e < nb_element; ++e) {
	Int * tmp_local_eq_nb_val = local_eq_nb_val;
	for (UInt i = 0; i < nb_nodes_per_element; ++i) {
	  UInt n = conn_val[i];
	  for (UInt d = 0; d < nb_degree_of_freedom; ++d) {
	    /**
	     * !!!!!!Careful!!!!!! This is a ugly fix. @todo this is a
	     * very ugly fix, because the offset for the global
	     * equation number, where the dof will be assembled, is
	     * hardcoded. In the future a class dof manager has to be
	     * added to Akantu to handle the mapping between the dofs
	     * and the equation numbers
	     * 
	     */

	    *tmp_local_eq_nb_val++ = eq_nb_val[n * nb_degree_of_freedom + d] - (dof_synchronizer.isPureGhostDOF(n * nb_degree_of_freedom + d) ? nb_global_dofs : 0);
	    
	  }
	}
	

	for (UInt i = 0; i < size_mat; ++i) {
	  Int c_irn = local_eq_nb_val[i];
	  UInt j_start = 0;
	  for (UInt j = j_start; j < size_mat; ++j) {
	    Int c_jcn = local_eq_nb_val[j];
	    Array<Int> index_pair(2,1);
	    index_pair(0) = c_irn;
	    index_pair(1) = c_jcn;
	    AOApplicationToPetsc(this->petsc_matrix_wrapper->ao, 2, index_pair.storage());
	    if (index_pair(0) >= first_global_index && index_pair(0) < first_global_index + this->local_size)  {
	      KeyCOO irn_jcn = keyPETSc(c_irn, c_jcn);
	      irn_jcn_k_it = irn_jcn_k.find(irn_jcn);

	      if (irn_jcn_k_it == irn_jcn_k.end()) {
		irn_jcn_k[irn_jcn] = nb_non_zero;
		
		
		// check if node is slave node
		if (index_pair(1) >= first_global_index && index_pair(1) < first_global_index + this->local_size)
		  this->d_nnz(index_pair(0) - first_global_index) += 1;
		else
		  this->o_nnz(index_pair(0) - first_global_index) += 1;
		nb_non_zero++;
		
	      }
	    }
	    
	  }
	  
	}
	conn_val += nb_nodes_per_element;
      }

      delete [] local_eq_nb_val;
    }
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
  


  // std::string mat_type;
  // mat_type.resize(20, 'a');
  // std::cout << "MatType: " << mat_type << std::endl;
  // const char * mat_type_ptr = mat_type.c_str();
  MatType type;
  MatGetType(this->petsc_matrix_wrapper->mat, &type);
  /// std::cout << "the matrix type is: " << type << std::endl;
  /**
   * PETSc will use only one of the following preallocation commands
   * depending on the matrix type and ignore the rest. Note that for
   * the block matrix format a block size of 1 is used. This might
   * result in bad performance. @todo For better results implement
   * buildProfile() with larger block size.
   * 
   */
  /// build profile:
  if (strcmp(type, MATSEQAIJ) == 0) {
    ierr = MatSeqAIJSetPreallocation(this->petsc_matrix_wrapper->mat,
				       0, d_nnz.storage()); CHKERRXX(ierr);
  } else if ((strcmp(type, MATMPIAIJ) == 0)) {
    ierr = MatMPIAIJSetPreallocation(this->petsc_matrix_wrapper->mat,
                                     0, d_nnz.storage(), 0,
                                     o_nnz.storage()); CHKERRXX(ierr);
  } else {
    AKANTU_DEBUG_ERROR("The type " << type << " of PETSc matrix is not handled by"
                       << " akantu in the preallocation step");
  }

  //ierr =  MatSeqSBAIJSetPreallocation(this->petsc_matrix_wrapper->mat, 1,
  //				      0, d_nnz.storage()); CHKERRXX(ierr);

  if (this->sparse_matrix_type==_symmetric) {
    /// set flag for symmetry to enable ICC/Cholesky preconditioner
    ierr = MatSetOption(this->petsc_matrix_wrapper->mat, MAT_SYMMETRIC, PETSC_TRUE); CHKERRXX(ierr);
    /// set flag for symmetric positive definite
    ierr = MatSetOption(this->petsc_matrix_wrapper->mat, MAT_SPD, PETSC_TRUE); CHKERRXX(ierr);
  }
  /// once the profile has been build ignore any new nonzero locations
  ierr = MatSetOption(this->petsc_matrix_wrapper->mat, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE); CHKERRXX(ierr);

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
/**
 * Method to save the nonzero pattern and the values stored at each position
 * @param filename name of the file in which the information will be stored
 */
void PETScMatrix::saveMatrix(const std::string & filename) const{
  AKANTU_DEBUG_IN();

  PetscErrorCode ierr;

  /// create Petsc viewer
  PetscViewer viewer; 
  ierr = PetscViewerASCIIOpen(this->petsc_matrix_wrapper->communicator, filename.c_str(), &viewer); CHKERRXX(ierr);

  /// set the format
  PetscViewerSetFormat(viewer, PETSC_VIEWER_DEFAULT); CHKERRXX(ierr);
  /// save the matrix
  /// @todo Write should be done in serial -> might cause problems
  ierr =  MatView(this->petsc_matrix_wrapper->mat, viewer); CHKERRXX(ierr);

  /// destroy the viewer
  ierr =  PetscViewerDestroy(&viewer); CHKERRXX(ierr);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Method to add an Akantu sparse matrix to the PETSc matrix
 * @param matrix Akantu sparse matrix to be added
 * @param alpha the factor specifying how many times the matrix should be added
 */
void PETScMatrix::add(const SparseMatrix & matrix, Real alpha) {
  PetscErrorCode ierr;
  //  AKANTU_DEBUG_ASSERT(nb_non_zero == matrix.getNbNonZero(),
  //		      "The two matrix don't have the same profiles");

  Real val_to_add = 0;
  for (UInt n = 0; n < matrix.getNbNonZero(); ++n) {
    Array<Int> index(2,1);
    UInt mat_to_add_offset = matrix.getOffset(); 
    index(0) = matrix.getIRN()(n)-mat_to_add_offset;
    index(1) =  matrix.getJCN()(n)-mat_to_add_offset;
    AOApplicationToPetsc(this->petsc_matrix_wrapper->ao, 2, index.storage());
    if (this->sparse_matrix_type == _symmetric && index(0) > index(1))
      std::swap(index(0), index(1));
    
    val_to_add = matrix.getA()(n) * alpha;
    /// MatSetValue might be very slow for MATBAIJ, might need to use MatSetValuesBlocked
    ierr = MatSetValue(this->petsc_matrix_wrapper->mat, index(0), index(1), val_to_add, ADD_VALUES); CHKERRXX(ierr);
    /// chek if sparse matrix to be added is symmetric. In this case
    /// the value also needs to be added at the transposed location in
    /// the matrix because PETSc is using the full profile, also for symmetric matrices
    if (matrix.getSparseMatrixType() == _symmetric && index(0) != index(1))
      ierr = MatSetValue(this->petsc_matrix_wrapper->mat, index(1), index(0), val_to_add, ADD_VALUES); CHKERRXX(ierr);
  }

  this->performAssembly();
}

/* -------------------------------------------------------------------------- */
/**
 * Method to add another PETSc matrix to this PETSc matrix
 * @param matrix PETSc matrix to be added
 * @param alpha the factor specifying how many times the matrix should be added
 */
void PETScMatrix::add(const PETScMatrix & matrix, Real alpha) {
  PetscErrorCode ierr;

  ierr = MatAXPY(this->petsc_matrix_wrapper->mat, alpha, matrix.petsc_matrix_wrapper->mat, SAME_NONZERO_PATTERN); CHKERRXX(ierr);


  this->performAssembly();
}

/* -------------------------------------------------------------------------- */
/**
 * MatSetValues() generally caches the values. The matrix is ready to
 * use only after MatAssemblyBegin() and MatAssemblyEnd() have been
 * called. (http://www.mcs.anl.gov/petsc/)
 */
void PETScMatrix::performAssembly() {
  PetscErrorCode ierr;
  ierr = MatAssemblyBegin(this->petsc_matrix_wrapper->mat, MAT_FINAL_ASSEMBLY); CHKERRXX(ierr);
  ierr = MatAssemblyEnd(this->petsc_matrix_wrapper->mat, MAT_FINAL_ASSEMBLY); CHKERRXX(ierr);
}

/* -------------------------------------------------------------------------- */
/**
 * Method is called when the static solver is destroyed, just before
 * PETSc is finalized. So all PETSc objects need to be destroyed at
 * this point.
 */
void PETScMatrix::beforeStaticSolverDestroy() {
  AKANTU_DEBUG_IN();

  try{
    this->destroyInternalData();
  } catch(...) {}

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
 /// destroy all the PETSc data structures used for this matrix
void PETScMatrix::destroyInternalData() {
  AKANTU_DEBUG_IN();

  if(this->is_petsc_matrix_initialized) {
  
    PetscErrorCode ierr;

    ierr = MatDestroy(&(this->petsc_matrix_wrapper->mat)); CHKERRXX(ierr);
    delete petsc_matrix_wrapper;

    this->is_petsc_matrix_initialized = false;
  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
/// access K(i, j). Works only for dofs on this processor!!!! 
Real PETScMatrix::operator()(UInt i, UInt j) const{
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(this->dof_synchronizer->isLocalOrMasterDOF(i) && this->dof_synchronizer->isLocalOrMasterDOF(j), "Operator works only for dofs on this processor");

  Array<Int> index(2,1);
  index(0) = this->dof_synchronizer->getDOFGlobalID(i);
  index(1) = this->dof_synchronizer->getDOFGlobalID(j);
  AOApplicationToPetsc(this->petsc_matrix_wrapper->ao, 2, index.storage());

  Real value = 0;

  PetscErrorCode ierr;
  /// @todo MatGetValue might be very slow for MATBAIJ, might need to use MatGetValuesBlocked
  ierr = MatGetValues(this->petsc_matrix_wrapper->mat, 1, &index(0), 1, &index(1), &value); CHKERRXX(ierr);


  AKANTU_DEBUG_OUT();


  return value;

}

/* -------------------------------------------------------------------------- */
/**
 * Apply Dirichlet boundary conditions by zeroing the rows and columns which correspond to blocked dofs
 * @param boundary array of booleans which are true if the dof at this position is blocked
 * @param block_val the value in the diagonal entry of blocked rows
 */
void PETScMatrix::applyBoundary(const Array<bool> & boundary, Real block_val) {
  AKANTU_DEBUG_IN();

  PetscErrorCode ierr;

  /// get the global equation numbers to find the rows that need to be zeroed for the blocked dofs
  Int * eq_nb_val = dof_synchronizer->getGlobalDOFEquationNumbers().storage();

  /// every processor calls the MatSetZero() only for his local or master dofs. This assures that not two processors or more try to zero the same row
  UInt nb_component = boundary.getNbComponent();
  UInt size = boundary.getSize();
  Int nb_blocked_local_master_eq_nb = 0;
  Array<Int> blocked_local_master_eq_nb(this->local_size);
  Int * blocked_lm_eq_nb_ptr  = blocked_local_master_eq_nb.storage();

  for (UInt i = 0; i < size; ++i) {
    for (UInt j = 0; j < nb_component; ++j) {
      UInt local_dof = i * nb_component + j;
      if (boundary(i, j) == true && this->dof_synchronizer->isLocalOrMasterDOF(local_dof)) {
	Int global_eq_nb = *eq_nb_val;
	*blocked_lm_eq_nb_ptr = global_eq_nb;
	++nb_blocked_local_master_eq_nb;
	++blocked_lm_eq_nb_ptr;
      }
      ++eq_nb_val;
    }
  }
  blocked_local_master_eq_nb.resize(nb_blocked_local_master_eq_nb);


  ierr = AOApplicationToPetsc(this->petsc_matrix_wrapper->ao, nb_blocked_local_master_eq_nb, blocked_local_master_eq_nb.storage() ); CHKERRXX(ierr);
  ierr = MatZeroRowsColumns(this->petsc_matrix_wrapper->mat, nb_blocked_local_master_eq_nb, blocked_local_master_eq_nb.storage(), block_val, 0, 0); CHKERRXX(ierr);

  this->performAssembly();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/// set all entries to zero while keeping the same nonzero pattern
void PETScMatrix::clear() {
  MatZeroEntries(this->petsc_matrix_wrapper->mat);
}

__END_AKANTU__
