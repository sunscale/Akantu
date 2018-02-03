/**
 * @file   sparse_matrix_petsc.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 * @date creation: Mon Dec 13 2010
 * @date last modification: Sat Feb 03 2018
 *
 * @brief  Implementation of PETSc matrix class
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "sparse_matrix_petsc.hh"
#include "communicator.hh"
#include "dof_manager_petsc.hh"
#include "mpi_communicator_data.hh"
/* -------------------------------------------------------------------------- */
#include <petscmat.h>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
SparseMatrixPETSc::SparseMatrixPETSc(DOFManagerPETSc & dof_manager,
                                     const MatrixType & sparse_matrix_type,
                                     const ID & id, const MemoryID & memory_id)
    : SparseMatrix(dof_manager, matrix_type, id, memory_id),
      dof_manager(dof_manager), {
  AKANTU_DEBUG_IN();

  auto mpi_comm = dof_manager.getMPIComm();

  PETSc_call(MatCreate, mpi_comm, &mat);

  if (sparse_matrix_type == _symmetric)
    PETSc_call(MatSetOption, MAT_SYMMETRIC, PETSC_TRUE);
  PETSc_call(MatSetOption, MAT_ROW_ORIENTED, PETSC_TRUE);

  PETSc_call(MatSetSizes, PETSC_DECIDE, PETSC_DECIDE, size_, size_);

  auto & is_ltog_mapping = dof_manager.getISLocalToGlobalMapping();

  PETSc_call(MatSetLocalToGlobalMapping, mat, is_ltog_mapping, is_ltog_mapping);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
SparseMatrixPETSc::~SparseMatrixPETSc() {
  AKANTU_DEBUG_IN();

  if (mat)
    PETSc_call(MatDestroy, &mat);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Method to save the nonzero pattern and the values stored at each position
 * @param filename name of the file in which the information will be stored
 */
void SparseMatrixPETSc::saveMatrix(const std::string & filename) const {
  AKANTU_DEBUG_IN();

  auto mpi_comm = dof_manager.getMPIComm();

  /// create Petsc viewer
  PetscViewer viewer;
  PETSc_call(PetscViewerASCIIOpen, mpi_comm, filename.c_str(), &viewer);

  /// set the format
  PETSc_call(PetscViewerSetFormat, viewer, PETSC_VIEWER_DEFAULT);

  /// save the matrix
  /// @todo Write should be done in serial -> might cause problems
  PETSc_call(MatView, mat, viewer);

  /// destroy the viewer
  PETSc_call(PetscViewerDestroy, &viewer);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/// Equivalent of *gemv in blas
void SparseMatrixPETSc::matVecMul(const Vec & x, Vec & y, Real alpha = 1.,
                                  Real beta = 0.) const {
  PETSc_call(MatMult, mat, x, y);

  if (alpha != 1.) {
    PETSc_call(VecScale, y, alpha);
  }

  PETSc_call(VecAYPX, y, beta, x);
}

/* -------------------------------------------------------------------------- */
void SparseMatrixPETSc::addMeToImpl(SparseMatrixPETSc & B, Real alpha) const {
  PETSc_call(MatAXPY, B.mat, alpha, mat, SAME_NONZERO_PATTERN);
  this->B.value_release++;
}

/* -------------------------------------------------------------------------- */
/**
 * Method to add another PETSc matrix to this PETSc matrix
 * @param matrix PETSc matrix to be added
 * @param alpha the factor specifying how many times the matrix should be added
 */
void SparseMatrixPETSc::addMeTo(SparseMatrix & B, Real alpha) const {
  if (auto * B_petsc = dynamic_cast<SparseMatrixPETSc *>(&B)) {
    this->addMeToImpl<SparseMatrixPETSc>(*B_petsc, alpha);
  } else {
    AKANTU_DEBUG_TO_IMPLEMENT();
    //    this->addMeToTemplated<SparseMatrix>(*this, alpha);
  }
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
  PETSc_call(MatAssemblyBegin, mat, MAT_FINAL_ASSEMBLY);
}

/* -------------------------------------------------------------------------- */
void SparseMatrixPETSc::endAssembly() {
  PETSc_call(MatAssemblyEnd, mat, MAT_FINAL_ASSEMBLY);
  PETSc_call(MatSetOption, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
  PETSc_call(MatSetOption, MAT_NEW_NONZERO_ALLOCATIONS, PETSC_TRUE);
  PETSc_call(MatSetOption, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
}

/* -------------------------------------------------------------------------- */
void SparseMatrixPETSc::applyBoundary(Real block_val) {
  AKANTU_DEBUG_IN();

  const auto & blocked_dofs = this->dof_manager.getGlobalBlockedDOFs();

  for (auto && data : arange(blocked_dofs)) {
    auto & dof = std::get<0>(data);
    auto & blocked = std::get<1>(data);

    if (blocked)
      blocked_list.push_back(dof);
  }

  PETSc_call(MatZeroRowsColumnsLocal, mat, blocked_list.size(),
             blocked_list.storage(), block_val, 0, 0);

  this->value_release++;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SparseMatrixPETSc::mul(Real alpha) {
  PETSc_call(MatScale, mat, alpha);
  this->value_release++;
}

/* -------------------------------------------------------------------------- */
void SparseMatrixPETSc::clear() {
  PETSc_call(MatZeroEntries, mat);
  this->value_release++;
}

/* -------------------------------------------------------------------------- */
void SparseMatrixPETSc::clearProfile() {
  SparseMatrix::clearProfile();
  PETSc_call(MatSetOption, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
  PETSc_call(MatSetOption, MAT_NEW_NONZERO_ALLOCATIONS, PETSC_TRUE);
  PETSc_call(MatSetOption, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);

  PETSc_call(MatZeroEntries, mat);
  this->value_release++;
}

/* -------------------------------------------------------------------------- */
void SparseMatrixPETSc::add(UInt i, UInt j) {
  PETSc_call(MatSetValue, i, j, 0, ADD_VALUES);
  this->value_release++;
}

/* -------------------------------------------------------------------------- */
void SparseMatrixPETSc::add(UInt i, UInt j, Real val) {
  PETSc_call(MatSetValue, i, j, val, ADD_VALUES);
  this->value_release++;
}

/* -------------------------------------------------------------------------- */
void SparseMatrixPETSc::addLocal(UInt i, UInt j) {
  PETSc_call(MatSetValueLocal, i, j, 0, ADD_VALUES);
  this->value_release++;
}

/* -------------------------------------------------------------------------- */
void SparseMatrixPETSc::addLocal(UInt i, UInt j, Real val) {
  PETSc_call(MatSetValueLocal, i, j, val, ADD_VALUES);
  this->value_release++;
}

/* -------------------------------------------------------------------------- */
void SparseMatrixPETSc::addLocal(const Vector<UInt> & rows,
                                 const Vector<UInt> & cols,
                                 const Matrix<Real> & vals) {
  PETSc_call(MatSetValuesLocal, rows.size(), rows.storage(), cols.size(),
             cols.storage(), vals.storage(), ADD_VALUES);
  this->value_release++;
}

/* -------------------------------------------------------------------------- */

} // akantu
