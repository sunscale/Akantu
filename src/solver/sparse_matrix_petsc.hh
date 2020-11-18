/**
 * @file   sparse_matrix_petsc.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Tue Feb 06 2018
 *
 * @brief  Interface for PETSc matrices
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

#ifndef AKANTU_PETSC_MATRIX_HH_
#define AKANTU_PETSC_MATRIX_HH_

/* -------------------------------------------------------------------------- */
#include "sparse_matrix.hh"
/* -------------------------------------------------------------------------- */
#include <petscmat.h>
/* -------------------------------------------------------------------------- */

namespace akantu {
class DOFManagerPETSc;
}

namespace akantu {

class SparseMatrixPETSc : public SparseMatrix {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  SparseMatrixPETSc(DOFManagerPETSc & dof_manager,
                    const MatrixType & matrix_type,
                    const ID & id = "sparse_matrix_petsc");

  SparseMatrixPETSc(const SparseMatrixPETSc & matrix,
                    const ID & id = "sparse_matrix_petsc");

  ~SparseMatrixPETSc() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// set the matrix to 0
  void zero() override;
  void set(Real /*val*/) override {
    AKANTU_TO_IMPLEMENT();
  }
  void clearProfile() override;

  /// add a non-zero element to the profile
  UInt add(UInt i, UInt j) override;

  /// assemble a local matrix in the sparse one
  void add(UInt i, UInt j, Real value) override;

  void addLocal(UInt i, UInt j);
  void addLocal(UInt i, UInt j, Real val);

  void addLocal(const Vector<Int> & rows, const Vector<Int> & cols,
                const Matrix<Real> & values);

  /// add a block of values
  void addValues(const Vector<Int> & rows, const Vector<Int> & cols,
                 const Matrix<Real> & values, MatrixType values_type);

  /// save the profil in a file using the MatrixMarket file format
  // void saveProfile(__attribute__((unused)) const std::string &) const
  // override {
  //   AKANTU_DEBUG_TO_IMPLEMENT();
  // }

  /// save the matrix in a file using the MatrixMarket file format
  void saveMatrix(const std::string & filename) const override;

  /// multiply the matrix by a coefficient
  void mul(Real alpha) override;

  /// Equivalent of *gemv in blas
  void matVecMul(const SolverVector & x, SolverVector & y, Real alpha = 1.,
                 Real beta = 0.) const override;

  /// modify the matrix to "remove" the blocked dof
  void applyBoundary(Real block_val = 1.) override;

  /// copy the profile of a matrix
  void copyProfile(const SparseMatrix & other) override;

  void applyModifications();

  void resize();

protected:
  void addMeTo(SparseMatrix & B, Real alpha) const override;

  /// This is the specific implementation
  void addMeToImpl(SparseMatrixPETSc & B, Real alpha) const;

  void beginAssembly();
  void endAssembly();

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// return the values at potition i, j
  inline Real operator()(UInt /*i*/, UInt /*j*/) const override {
    AKANTU_TO_IMPLEMENT();
  }
  /// return the values at potition i, j
  inline Real & operator()(UInt /*i*/, UInt /*j*/) override {
    AKANTU_TO_IMPLEMENT();
  }

  UInt getRelease() const override { return release; };

  operator Mat &() { return mat; }
  operator const Mat &() const { return mat; }
  AKANTU_GET_MACRO(Mat, mat, const Mat &);
  AKANTU_GET_MACRO_NOT_CONST(Mat, mat, Mat &);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  // DOFManagerPETSc that contains the numbering for petsc
  DOFManagerPETSc & dof_manager;

  /// store the PETSc matrix
  Mat mat;

  /// matrix release
  UInt release{0};
};

} // namespace akantu

#endif /* AKANTU_PETSC_MATRIX_HH_ */
