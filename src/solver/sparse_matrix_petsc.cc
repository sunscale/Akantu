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
#include "sparse_matrix_petsc.hh"
#include "mpi_type_wrapper.hh"
#include "dof_manager_petsc.hh"
#include "static_communicator.hh"
/* -------------------------------------------------------------------------- */
#include <cstring>
#include <petscsys.h>
/* -------------------------------------------------------------------------- */

namespace akantu {

#if not defined(PETSC_CLANGUAGE_CXX)
int aka_PETScError(int ierr) {
  CHKERRQ(ierr);
  return 0;
}
#endif

/* -------------------------------------------------------------------------- */
SparseMatrixPETSc::SparseMatrixPETSc(DOFManagerPETSc & dof_manager,
    const MatrixType & sparse_matrix_type, const ID & id,
    const MemoryID & memory_id)
  : SparseMatrix(dof_manager, matrix_type, id, memory_id),
    dof_manager(dof_manager),
    d_nnz(0, 1, "dnnz"),
    o_nnz(0, 1, "onnz"), first_global_index(0) {
  AKANTU_DEBUG_IN();

  PetscErrorCode ierr;

  // create the PETSc matrix object
  ierr = MatCreate(PETSC_COMM_WORLD, &this->mat);
  CHKERRXX(ierr);

  /**
   * Set the matrix type
   * @todo PETSc does currently not support a straightforward way to
   * apply Dirichlet boundary conditions for MPISBAIJ
   * matrices. Therefore always the entire matrix is allocated. It
   * would be possible to use MATSBAIJ for sequential matrices in case
   * memory becomes critical. Also, block matrices would give a much
   * better performance. Modify this in the future!
   */
  ierr = MatSetType(this->mat, MATAIJ);
  CHKERRXX(ierr);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
SparseMatrixPETSc::~SparseMatrixPETSc() {
  AKANTU_DEBUG_IN();

  /// destroy all the PETSc data structures used for this matrix
  PetscErrorCode ierr;
  ierr = MatDestroy(&this->mat);
  CHKERRXX(ierr);

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
void SparseMatrixPETSc::setSize() {
  AKANTU_DEBUG_IN();

  //  PetscErrorCode ierr;

  /// find the number of dofs corresponding to master or local nodes
  UInt nb_dofs = this->dof_manager.getLocalSystemSize();
  // UInt nb_local_master_dofs = 0;

  /// create array to store the global equation number of all local and master
  /// dofs
  Array<Int> local_master_eq_nbs(nb_dofs);
  Array<Int>::scalar_iterator it_eq_nb = local_master_eq_nbs.begin();

  throw;
  /// get the pointer to the global equation number array
 //  Int * eq_nb_val =
//       this->dof_synchronizer->getGlobalDOFEquationNumbers().storage();

//   for (UInt i = 0; i < nb_dofs; ++i) {
//     if (this->dof_synchronizer->isLocalOrMasterDOF(i)) {
//       *it_eq_nb = eq_nb_val[i];
//       ++it_eq_nb;
//       ++nb_local_master_dofs;
//     }
//   }

//   local_master_eq_nbs.resize(nb_local_master_dofs);

//   /// set the local size
//   this->local_size = nb_local_master_dofs;

// /// resize PETSc matrix
// #if defined(AKANTU_USE_MPI)
//   ierr = MatSetSizes(this->petsc_matrix_wrapper->mat, this->local_size,
//                      this->local_size, this->size, this->size);
//   CHKERRXX(ierr);
// #else
//   ierr = MatSetSizes(this->petsc_matrix_wrapper->mat, this->local_size,
//                      this->local_size);
//   CHKERRXX(ierr);
// #endif

//   /// create mapping from akantu global numbering to petsc global numbering
//   this->createGlobalAkantuToPETScMap(local_master_eq_nbs.storage());

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
void
SparseMatrixPETSc::createGlobalAkantuToPETScMap(Int * local_master_eq_nbs_ptr) {
  AKANTU_DEBUG_IN();

  PetscErrorCode ierr;

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  UInt rank = comm.whoAmI();

  // initialize vector to store the number of local and master nodes that are
  // assigned to each processor
  Vector<UInt> master_local_ndofs_per_proc(nb_proc);

  /// store the nb of master and local dofs on each processor
  master_local_ndofs_per_proc(rank) = this->local_size;

  /// exchange the information among all processors
  comm.allGather(master_local_ndofs_per_proc.storage(), 1);

  /// each processor creates a map for his akantu global dofs to the
  /// corresponding petsc global dofs

  /// determine the PETSc-index for the first dof on each processor

  for (UInt i = 0; i < rank; ++i) {
    this->first_global_index += master_local_ndofs_per_proc(i);
  }

  /// create array for petsc ordering
  Array<Int> petsc_dofs(this->local_size);
  Array<Int>::scalar_iterator it_petsc_dofs = petsc_dofs.begin();

  for (Int i = this->first_global_index;
       i < this->first_global_index + this->local_size; ++i, ++it_petsc_dofs) {
    *it_petsc_dofs = i;
  }

  ierr = AOCreateBasic(PETSC_COMM_WORLD,
                       this->local_size, local_master_eq_nbs_ptr,
                       petsc_dofs.storage(), &(this->ao));
  CHKERRXX(ierr);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Method to save the nonzero pattern and the values stored at each position
 * @param filename name of the file in which the information will be stored
 */
void SparseMatrixPETSc::saveMatrix(const std::string & filename) const {
  AKANTU_DEBUG_IN();

  PetscErrorCode ierr;

  /// create Petsc viewer
  PetscViewer viewer;
  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,
                              filename.c_str(), &viewer);
  CHKERRXX(ierr);

  /// set the format
  PetscViewerSetFormat(viewer, PETSC_VIEWER_DEFAULT);
  CHKERRXX(ierr);

  /// save the matrix
  /// @todo Write should be done in serial -> might cause problems
  ierr = MatView(this->mat, viewer);
  CHKERRXX(ierr);

  /// destroy the viewer
  ierr = PetscViewerDestroy(&viewer);
  CHKERRXX(ierr);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Method to add an Akantu sparse matrix to the PETSc matrix
 * @param matrix Akantu sparse matrix to be added
 * @param alpha the factor specifying how many times the matrix should be added
 */
// void SparseMatrixPETSc::add(const SparseMatrix & matrix, Real alpha) {
//   PetscErrorCode ierr;
//   //  AKANTU_DEBUG_ASSERT(nb_non_zero == matrix.getNbNonZero(),
//   //                  "The two matrix don't have the same profiles");

//   Real val_to_add = 0;
//   Array<Int> index(2);
//   for (UInt n = 0; n < matrix.getNbNonZero(); ++n) {
//     UInt mat_to_add_offset = matrix.getOffset();
//     index(0) = matrix.getIRN()(n) - mat_to_add_offset;
//     index(1) = matrix.getJCN()(n) - mat_to_add_offset;
//     AOApplicationToPetsc(this->petsc_matrix_wrapper->ao, 2, index.storage());
//     if (this->sparse_matrix_type == _symmetric && index(0) > index(1))
//       std::swap(index(0), index(1));

//     val_to_add = matrix.getA()(n) * alpha;
//     /// MatSetValue might be very slow for MATBAIJ, might need to use
//     /// MatSetValuesBlocked
//     ierr = MatSetValue(this->petsc_matrix_wrapper->mat, index(0), index(1),
//                        val_to_add, ADD_VALUES);
//     CHKERRXX(ierr);
//     /// chek if sparse matrix to be added is symmetric. In this case
//     /// the value also needs to be added at the transposed location in
//     /// the matrix because PETSc is using the full profile, also for symmetric
//     /// matrices
//     if (matrix.getSparseMatrixType() == _symmetric && index(0) != index(1))
//       ierr = MatSetValue(this->petsc_matrix_wrapper->mat, index(1), index(0),
//                          val_to_add, ADD_VALUES);
//     CHKERRXX(ierr);
//   }

//   this->performAssembly();
// }

/* -------------------------------------------------------------------------- */
/**
 * Method to add another PETSc matrix to this PETSc matrix
 * @param matrix PETSc matrix to be added
 * @param alpha the factor specifying how many times the matrix should be added
 */
void SparseMatrixPETSc::add(const SparseMatrixPETSc & matrix, Real alpha) {
  PetscErrorCode ierr;

  ierr = MatAXPY(this->mat, alpha,
                 matrix.mat, SAME_NONZERO_PATTERN);
  CHKERRXX(ierr);

  this->performAssembly();
}

/* -------------------------------------------------------------------------- */
/**
 * MatSetValues() generally caches the values. The matrix is ready to
 * use only after MatAssemblyBegin() and MatAssemblyEnd() have been
 * called. (http://www.mcs.anl.gov/petsc/)
 */
void SparseMatrixPETSc::performAssembly() {
  this->beginAssembly();
  this->endAssembly();
}

/* -------------------------------------------------------------------------- */
void SparseMatrixPETSc::beginAssembly() {
  PetscErrorCode ierr;
  ierr = MatAssemblyBegin(this->mat, MAT_FINAL_ASSEMBLY);
  CHKERRXX(ierr);
  ierr = MatAssemblyEnd(this->mat, MAT_FINAL_ASSEMBLY);
  CHKERRXX(ierr);
}

/* -------------------------------------------------------------------------- */
void SparseMatrixPETSc::endAssembly() {
  PetscErrorCode ierr;
  ierr = MatAssemblyEnd(this->mat, MAT_FINAL_ASSEMBLY);
  CHKERRXX(ierr);
}

/* -------------------------------------------------------------------------- */
/// access K(i, j). Works only for dofs on this processor!!!!
Real SparseMatrixPETSc::operator()(UInt i, UInt j) const {
  AKANTU_DEBUG_IN();

  // AKANTU_DEBUG_ASSERT(this->dof_synchronizer->isLocalOrMasterDOF(i) &&
  //                         this->dof_synchronizer->isLocalOrMasterDOF(j),
  //                     "Operator works only for dofs on this processor");

  // Array<Int> index(2, 1);
  // index(0) = this->dof_synchronizer->getDOFGlobalID(i);
  // index(1) = this->dof_synchronizer->getDOFGlobalID(j);
  // AOApplicationToPetsc(this->petsc_matrix_wrapper->ao, 2, index.storage());

  // Real value = 0;

  // PetscErrorCode ierr;
  // /// @todo MatGetValue might be very slow for MATBAIJ, might need to use
  // /// MatGetValuesBlocked
  // ierr = MatGetValues(this->petsc_matrix_wrapper->mat, 1, &index(0), 1,
  //                     &index(1), &value);
  // CHKERRXX(ierr);

  // AKANTU_DEBUG_OUT();

  // return value;
  return 0.;
}

/* -------------------------------------------------------------------------- */
/**
 * Apply Dirichlet boundary conditions by zeroing the rows and columns which
 * correspond to blocked dofs
 * @param boundary array of booleans which are true if the dof at this position
 * is blocked
 * @param block_val the value in the diagonal entry of blocked rows
 */
void SparseMatrixPETSc::applyBoundary(const Array<bool> & boundary,
                                      Real block_val) {
  AKANTU_DEBUG_IN();

  // PetscErrorCode ierr;

  // /// get the global equation numbers to find the rows that need to be zeroed
  // /// for the blocked dofs
  // Int * eq_nb_val = dof_synchronizer->getGlobalDOFEquationNumbers().storage();

  // /// every processor calls the MatSetZero() only for his local or master dofs.
  // /// This assures that not two processors or more try to zero the same row
  // UInt nb_component = boundary.getNbComponent();
  // UInt size = boundary.size();
  // Int nb_blocked_local_master_eq_nb = 0;
  // Array<Int> blocked_local_master_eq_nb(this->local_size);
  // Int * blocked_lm_eq_nb_ptr = blocked_local_master_eq_nb.storage();

  // for (UInt i = 0; i < size; ++i) {
  //   for (UInt j = 0; j < nb_component; ++j) {
  //     UInt local_dof = i * nb_component + j;
  //     if (boundary(i, j) == true &&
  //         this->dof_synchronizer->isLocalOrMasterDOF(local_dof)) {
  //       Int global_eq_nb = *eq_nb_val;
  //       *blocked_lm_eq_nb_ptr = global_eq_nb;
  //       ++nb_blocked_local_master_eq_nb;
  //       ++blocked_lm_eq_nb_ptr;
  //     }
  //     ++eq_nb_val;
  //   }
  // }
  // blocked_local_master_eq_nb.resize(nb_blocked_local_master_eq_nb);

  // ierr = AOApplicationToPetsc(this->petsc_matrix_wrapper->ao,
  //                             nb_blocked_local_master_eq_nb,
  //                             blocked_local_master_eq_nb.storage());
  // CHKERRXX(ierr);
  // ierr = MatZeroRowsColumns(
  //     this->petsc_matrix_wrapper->mat, nb_blocked_local_master_eq_nb,
  //     blocked_local_master_eq_nb.storage(), block_val, 0, 0);
  // CHKERRXX(ierr);

  // this->performAssembly();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/// set all entries to zero while keeping the same nonzero pattern
void SparseMatrixPETSc::clear() {
  MatZeroEntries(this->mat);
}

} // akantu
