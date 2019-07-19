/**
 * @file   sparse_matrix.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Dec 13 2010
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  sparse matrix storage class (distributed assembled matrix)
 * This is a COO format (Coordinate List)
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
#include "aka_common.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_SPARSE_MATRIX_HH__
#define __AKANTU_SPARSE_MATRIX_HH__

/* -------------------------------------------------------------------------- */
namespace akantu {
class DOFManager;
class TermsToAssemble;
class SolverVector;
} // namespace akantu

namespace akantu {

class SparseMatrix {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  SparseMatrix(DOFManager & dof_manager, const MatrixType & matrix_type,
               const ID & id = "sparse_matrix");

  SparseMatrix(const SparseMatrix & matrix, const ID & id = "sparse_matrix");

  virtual ~SparseMatrix();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// remove the existing profile
  virtual void clearProfile();

  /// set the matrix to 0
  virtual void clear() = 0;

  /// add a non-zero element to the profile
  virtual UInt add(UInt i, UInt j) = 0;

  /// assemble a local matrix in the sparse one
  virtual void add(UInt i, UInt j, Real value) = 0;

  /// save the profil in a file using the MatrixMarket file format
  virtual void saveProfile(__attribute__((unused)) const std::string &) const {
    AKANTU_TO_IMPLEMENT();
  }

  /// save the matrix in a file using the MatrixMarket file format
  virtual void saveMatrix(__attribute__((unused)) const std::string &) const {
    AKANTU_TO_IMPLEMENT();
  };

  /// multiply the matrix by a coefficient
  virtual void mul(Real alpha) = 0;

  /// add matrices
  virtual void add(const SparseMatrix & matrix, Real alpha = 1.);

  /// Equivalent of *gemv in blas
  virtual void matVecMul(const SolverVector & x, SolverVector & y,
                         Real alpha = 1., Real beta = 0.) const = 0;

  /// modify the matrix to "remove" the blocked dof
  virtual void applyBoundary(Real block_val = 1.) = 0;

  /// copy the profile of another matrix
  virtual void copyProfile(const SparseMatrix & other) = 0;

  /// operator *=
  SparseMatrix & operator*=(Real alpha) {
    this->mul(alpha);
    return *this;
  }

protected:
  /// This is the revert of add B += \alpha * *this;
  virtual void addMeTo(SparseMatrix & B, Real alpha) const = 0;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// return the values at potition i, j
  virtual inline Real operator()(__attribute__((unused)) UInt i,
                                 __attribute__((unused)) UInt j) const {
    AKANTU_TO_IMPLEMENT();
  }
  /// return the values at potition i, j
  virtual inline Real & operator()(__attribute__((unused)) UInt i,
                                   __attribute__((unused)) UInt j) {
    AKANTU_TO_IMPLEMENT();
  }

  AKANTU_GET_MACRO(NbNonZero, nb_non_zero, UInt);
  UInt size() const { return size_; }
  AKANTU_GET_MACRO(MatrixType, matrix_type, const MatrixType &);

  virtual UInt getRelease() const = 0;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  ID id;

  /// Underlying dof manager
  DOFManager & _dof_manager;

  /// sparce matrix type
  MatrixType matrix_type;

  /// Size of the matrix
  UInt size_;

  /// number of processors
  UInt nb_proc;

  /// number of non zero element
  UInt nb_non_zero;
};

// Array<Real> & operator*=(Array<Real> & vect, const SparseMatrix & mat);

} // namespace akantu

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "sparse_matrix_inline_impl.cc"

#endif /* __AKANTU_SPARSE_MATRIX_HH__ */
