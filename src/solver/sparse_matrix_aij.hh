/**
 * @file   sparse_matrix_aij.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Dec 13 2010
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  AIJ implementation of the SparseMatrix (this the format used by
 * Mumps)
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
#include "aka_array.hh"
#include "aka_common.hh"
#include "sparse_matrix.hh"
/* -------------------------------------------------------------------------- */
#include <unordered_map>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_SPARSE_MATRIX_AIJ_HH_
#define AKANTU_SPARSE_MATRIX_AIJ_HH_

namespace akantu {
class DOFManagerDefault;
class TermsToAssemble;
} // namespace akantu

namespace akantu {

class SparseMatrixAIJ : public SparseMatrix {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  SparseMatrixAIJ(DOFManagerDefault & dof_manager,
                  const MatrixType & matrix_type,
                  const ID & id = "sparse_matrix_aij");

  SparseMatrixAIJ(const SparseMatrixAIJ & matrix,
                  const ID & id = "sparse_matrix_aij");

  ~SparseMatrixAIJ() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// remove the existing profile
  inline void clearProfile() override;

  /// add a non-zero element
  inline UInt add(UInt i, UInt j) override;

  /// set the matrix to 0
  void set(Real val) override;

  /// assemble a local matrix in the sparse one
  inline void add(UInt i, UInt j, Real value) override;

  /// add a block of values
  inline void addValues(const Vector<Int> & is, const Vector<Int> & js,
                        const Matrix<Real> & values, MatrixType values_type);

  /// set the size of the matrix
  void resize(UInt size) { this->size_ = size; }

  /// modify the matrix to "remove" the blocked dof
  void applyBoundary(Real block_val = 1.) override;

  /// save the profil in a file using the MatrixMarket file format
  void saveProfile(const std::string & filename) const override;

  /// save the matrix in a file using the MatrixMarket file format
  void saveMatrix(const std::string & filename) const override;

  /// copy assuming the profile are the same
  virtual void copyContent(const SparseMatrix & matrix);

  /// multiply the matrix by a scalar
  void mul(Real alpha) override;

  /// Equivalent of *gemv in blas
  void matVecMul(const SolverVector & x, SolverVector & y, Real alpha = 1.,
                 Real beta = 0.) const override;

  void matVecMul(const Array<Real> & x, Array<Real> & y, Real alpha = 1.,
                 Real beta = 0.) const;

  /// copy the profile of another matrix
  void copyProfile(const SparseMatrix & other) override;

  /* ------------------------------------------------------------------------ */
  /// accessor to A_{ij} - if (i, j) not present it returns 0
  inline Real operator()(UInt i, UInt j) const override;

  /// accessor to A_{ij} - if (i, j) not present it fails, (i, j) should be
  /// first added to the profile
  inline Real & operator()(UInt i, UInt j) override;

protected:
  void addMeTo(SparseMatrix & B, Real alpha) const override;

  inline void addSymmetricValuesToSymmetric(const Vector<Int> & is,
                                            const Vector<Int> & js,
                                            const Matrix<Real> & values);
  inline void addUnsymmetricValuesToSymmetric(const Vector<Int> & is,
                                              const Vector<Int> & js,
                                              const Matrix<Real> & values);
  inline void addValuesToUnsymmetric(const Vector<Int> & is,
                                     const Vector<Int> & js,
                                     const Matrix<Real> & values);

private:
  /// This is just to inline the addToMatrix function
  template <class MatrixType>
  void addMeToTemplated(MatrixType & B, Real alpha) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(IRN, irn, const Array<Int> &);

  AKANTU_GET_MACRO(JCN, jcn, const Array<Int> &);

  AKANTU_GET_MACRO(A, a, const Array<Real> &);

  /// The release changes at each call of a function that changes the profile,
  /// it in increasing but could overflow so it should be checked as
  /// (my_release != release) and not as (my_release < release)
  AKANTU_GET_MACRO(ProfileRelease, profile_release, UInt);
  AKANTU_GET_MACRO(ValueRelease, value_release, UInt);
  UInt getRelease() const override { return value_release; }

protected:
  using KeyCOO = std::pair<UInt, UInt>;
  using coordinate_list_map = std::unordered_map<KeyCOO, UInt>;

  /// get the pair corresponding to (i, j)
  inline KeyCOO key(UInt i, UInt j) const {
    if (this->matrix_type == _symmetric && (i > j)) {
      return std::make_pair(j, i);
    }
    return std::make_pair(i, j);
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  DOFManagerDefault & dof_manager;

  /// row indexes
  Array<Int> irn;

  /// column indexes
  Array<Int> jcn;

  /// values : A[k] = Matrix[irn[k]][jcn[k]]
  Array<Real> a;

  /// Profile release
  UInt profile_release{1};

  /// Value release
  UInt value_release{1};

  /// map for (i, j) ->  k correspondence
  coordinate_list_map irn_jcn_k;
};

} // namespace akantu

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "sparse_matrix_aij_inline_impl.hh"

#endif /* AKANTU_SPARSE_MATRIX_AIJ_HH_ */
