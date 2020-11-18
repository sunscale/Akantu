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
#include "dof_manager_petsc.hh"
#include "mpi_communicator_data.hh"
#include "solver_vector_petsc.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
SparseMatrixPETSc::SparseMatrixPETSc(DOFManagerPETSc & dof_manager,
                                     const MatrixType & matrix_type,
                                     const ID & id)
    : SparseMatrix(dof_manager, matrix_type, id), dof_manager(dof_manager) {
  AKANTU_DEBUG_IN();

  auto && mpi_comm = dof_manager.getMPIComm();

  PETSc_call(MatCreate, mpi_comm, &mat);
  detail::PETScSetName(mat, id);

  resize();

  PETSc_call(MatSetFromOptions, mat);

  PETSc_call(MatSetUp, mat);

  PETSc_call(MatSetOption, mat, MAT_ROW_ORIENTED, PETSC_TRUE);
  PETSc_call(MatSetOption, mat, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);

  if (matrix_type == _symmetric) {
    PETSc_call(MatSetOption, mat, MAT_SYMMETRIC, PETSC_TRUE);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
SparseMatrixPETSc::SparseMatrixPETSc(const SparseMatrixPETSc & matrix,
                                     const ID & id)
    : SparseMatrix(matrix, id), dof_manager(matrix.dof_manager) {
  PETSc_call(MatDuplicate, matrix.mat, MAT_COPY_VALUES, &mat);
  detail::PETScSetName(mat, id);
}

/* -------------------------------------------------------------------------- */
SparseMatrixPETSc::~SparseMatrixPETSc() {
  AKANTU_DEBUG_IN();

  if (mat != nullptr) {
    PETSc_call(MatDestroy, &mat);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SparseMatrixPETSc::resize() {
  auto local_size = dof_manager.getPureLocalSystemSize();
  PETSc_call(MatSetSizes, mat, local_size, local_size, size_, size_);

  auto & is_ltog_mapping = dof_manager.getISLocalToGlobalMapping();
  PETSc_call(MatSetLocalToGlobalMapping, mat, is_ltog_mapping, is_ltog_mapping);
}

/* -------------------------------------------------------------------------- */
/**
 * Method to save the nonzero pattern and the values stored at each position
 * @param filename name of the file in which the information will be stored
 */
void SparseMatrixPETSc::saveMatrix(const std::string & filename) const {
  AKANTU_DEBUG_IN();

  auto && mpi_comm = dof_manager.getMPIComm();

  /// create Petsc viewer
  PetscViewer viewer;
  PETSc_call(PetscViewerASCIIOpen, mpi_comm, filename.c_str(), &viewer);
  PETSc_call(PetscViewerPushFormat, viewer, PETSC_VIEWER_ASCII_MATRIXMARKET);
  PETSc_call(MatView, mat, viewer);
  PETSc_call(PetscViewerPopFormat, viewer);
  PETSc_call(PetscViewerDestroy, &viewer);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/// Equivalent of *gemv in blas
void SparseMatrixPETSc::matVecMul(const SolverVector & _x, SolverVector & _y,
                                  Real alpha, Real beta) const {
  const auto & x = aka::as_type<SolverVectorPETSc>(_x);
  auto & y = aka::as_type<SolverVectorPETSc>(_y);

  // y = alpha A x + beta y
  SolverVectorPETSc w(x, this->id + ":tmp");

  // w = A x
  if (release == 0) {
    PETSc_call(VecZeroEntries, w);
  } else {
    PETSc_call(MatMult, mat, x, w);
  }

  if (alpha != 1.) {
    // w = alpha w
    PETSc_call(VecScale, w, alpha);
  }

  // y = w + beta y
  PETSc_call(VecAYPX, y, beta, w);
}

/* -------------------------------------------------------------------------- */
void SparseMatrixPETSc::addMeToImpl(SparseMatrixPETSc & B, Real alpha) const {
  PETSc_call(MatAXPY, B.mat, alpha, mat, SAME_NONZERO_PATTERN);

  B.release++;
}

/* -------------------------------------------------------------------------- */
void SparseMatrixPETSc::addMeTo(SparseMatrix & B, Real alpha) const {
  if (aka::is_of_type<SparseMatrixPETSc>(B)) {
    auto & B_petsc = aka::as_type<SparseMatrixPETSc>(B);
    this->addMeToImpl(B_petsc, alpha);
  } else {
    AKANTU_TO_IMPLEMENT();
    //    this->addMeToTemplated<SparseMatrix>(*this, alpha);
  }
}

/* -------------------------------------------------------------------------- */
/**
 * MatSetValues() generally caches the values. The matrix is ready to
 * use only after MatAssemblyBegin() and MatAssemblyEnd() have been
 * called. (http://www.mcs.anl.gov/petsc/)
 */
void SparseMatrixPETSc::applyModifications() {
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
  PETSc_call(MatSetOption, mat, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);

  this->release++;
}

/* -------------------------------------------------------------------------- */
void SparseMatrixPETSc::copyProfile(const SparseMatrix & other) {
  const auto & A = aka::as_type<SparseMatrixPETSc>(other);

  MatDestroy(&mat);
  MatDuplicate(A.mat, MAT_DO_NOT_COPY_VALUES, &mat);
}

/* -------------------------------------------------------------------------- */
void SparseMatrixPETSc::applyBoundary(Real block_val) {
  AKANTU_DEBUG_IN();

  const auto & blocked_dofs = this->dof_manager.getGlobalBlockedDOFs();
  // std::vector<PetscInt> rows;
  // for (auto && data : enumerate(blocked)) {
  //   if (std::get<1>(data)) {
  //     rows.push_back(std::get<0>(data));
  //   }
  // }
  // applyModifications();

  static int c = 0;

  saveMatrix("before_blocked_" + std::to_string(c) + ".mtx");

  PETSc_call(MatZeroRowsColumnsLocal, mat, blocked_dofs.size(),
             blocked_dofs.storage(), block_val, nullptr, nullptr);

  saveMatrix("after_blocked_" + std::to_string(c) + ".mtx");
  ++c;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SparseMatrixPETSc::mul(Real alpha) {
  PETSc_call(MatScale, mat, alpha);
  this->release++;
}

/* -------------------------------------------------------------------------- */
void SparseMatrixPETSc::zero() {
  PETSc_call(MatZeroEntries, mat);
  this->release++;
}

/* -------------------------------------------------------------------------- */
void SparseMatrixPETSc::clearProfile() {
  SparseMatrix::clearProfile();
  PETSc_call(MatResetPreallocation, mat);
  PETSc_call(MatSetOption, mat, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);
  //   PETSc_call(MatSetOption, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
  //   PETSc_call(MatSetOption, MAT_NEW_NONZERO_ALLOCATIONS, PETSC_TRUE);
  //   PETSc_call(MatSetOption, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);

  this->zero();
}

/* -------------------------------------------------------------------------- */
UInt SparseMatrixPETSc::add(UInt i, UInt j) {
  PETSc_call(MatSetValue, mat, i, j, 0, ADD_VALUES);
  return 0;
}

/* -------------------------------------------------------------------------- */
void SparseMatrixPETSc::add(UInt i, UInt j, Real val) {
  PETSc_call(MatSetValue, mat, i, j, val, ADD_VALUES);
}

/* -------------------------------------------------------------------------- */
void SparseMatrixPETSc::addLocal(UInt i, UInt j) {
  PETSc_call(MatSetValueLocal, mat, i, j, 0, ADD_VALUES);
}

/* -------------------------------------------------------------------------- */
void SparseMatrixPETSc::addLocal(UInt i, UInt j, Real val) {
  PETSc_call(MatSetValueLocal, mat, i, j, val, ADD_VALUES);
}

/* -------------------------------------------------------------------------- */
void SparseMatrixPETSc::addLocal(const Vector<Int> & rows,
                                 const Vector<Int> & cols,
                                 const Matrix<Real> & values) {
  PETSc_call(MatSetValuesLocal, mat, rows.size(), rows.storage(), cols.size(),
             cols.storage(), values.storage(), ADD_VALUES);
}

/* -------------------------------------------------------------------------- */
void SparseMatrixPETSc::addValues(const Vector<Int> & rows,
                                  const Vector<Int> & cols,
                                  const Matrix<Real> & values, MatrixType values_type) {
  if (values_type == _unsymmetric and matrix_type == _symmetric) {
    PETSc_call(MatSetOption, mat, MAT_SYMMETRIC, PETSC_FALSE);
    PETSc_call(MatSetOption, mat, MAT_STRUCTURALLY_SYMMETRIC, PETSC_FALSE);
  }

  PETSc_call(MatSetValues, mat, rows.size(), rows.storage(), cols.size(),
             cols.storage(), values.storage(), ADD_VALUES);
}

/* -------------------------------------------------------------------------- */

} // namespace akantu
