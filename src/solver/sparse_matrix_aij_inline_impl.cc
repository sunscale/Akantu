/**
 * @file   sparse_matrix_aij_inline_impl.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Aug 21 2015
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  Implementation of inline functions of SparseMatrixAIJ
 *
 * @section LICENSE
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "sparse_matrix_aij.hh"

/* -------------------------------------------------------------------------- */
#ifndef __AKANTU_SPARSE_MATRIX_AIJ_INLINE_IMPL_CC__
#define __AKANTU_SPARSE_MATRIX_AIJ_INLINE_IMPL_CC__

namespace akantu {

inline UInt SparseMatrixAIJ::add(UInt i, UInt j) {
  KeyCOO jcn_irn = this->key(i, j);

  auto it = this->irn_jcn_k.find(jcn_irn);

  if (!(it == this->irn_jcn_k.end()))
    return it->second;

  if (i + 1 > this->size_)
    this->size_ = i + 1;
  if (j + 1 > this->size_)
    this->size_ = j + 1;

  this->irn.push_back(i + 1);
  this->jcn.push_back(j + 1);
  this->a.push_back(0.);

  this->irn_jcn_k[jcn_irn] = this->nb_non_zero;

  (this->nb_non_zero)++;

  this->profile_release++;
  this->value_release++;

  return (this->nb_non_zero - 1);
}

/* -------------------------------------------------------------------------- */
inline void SparseMatrixAIJ::clearProfile() {
  SparseMatrix::clearProfile();

  this->irn_jcn_k.clear();

  this->irn.resize(0);
  this->jcn.resize(0);
  this->a.resize(0);

  this->size_ = 0;
  this->nb_non_zero = 0;

  this->profile_release++;
  this->value_release++;
}

/* -------------------------------------------------------------------------- */
inline void SparseMatrixAIJ::add(UInt i, UInt j, Real value) {
  UInt idx = this->add(i, j);

  this->a(idx) += value;

  this->value_release++;
}

/* -------------------------------------------------------------------------- */
inline Real SparseMatrixAIJ::operator()(UInt i, UInt j) const {
  KeyCOO jcn_irn = this->key(i, j);
  auto irn_jcn_k_it = this->irn_jcn_k.find(jcn_irn);

  if (irn_jcn_k_it == this->irn_jcn_k.end())
    return 0.;
  return this->a(irn_jcn_k_it->second);
}

/* -------------------------------------------------------------------------- */
inline Real & SparseMatrixAIJ::operator()(UInt i, UInt j) {
  KeyCOO jcn_irn = this->key(i, j);
  auto irn_jcn_k_it = this->irn_jcn_k.find(jcn_irn);
  AKANTU_DEBUG_ASSERT(irn_jcn_k_it != this->irn_jcn_k.end(),
                      "Couple (i,j) = (" << i << "," << j
                                         << ") does not exist in the profile");

  // it may change the profile so it is considered as a change
  this->value_release++;

  return this->a(irn_jcn_k_it->second);
}

/* -------------------------------------------------------------------------- */
inline void
SparseMatrixAIJ::addSymmetricValuesToSymmetric(const Vector<Int> & is,
                                               const Vector<Int> & js,
                                               const Matrix<Real> & values) {
  for (UInt i = 0; i < values.rows(); ++i) {
    UInt c_irn = is(i);
    if (c_irn < size_) {
      for (UInt j = i; j < values.cols(); ++j) {
        UInt c_jcn = js(j);
        if (c_jcn < size_) {
          operator()(c_irn, c_jcn) += values(i, j);
        }
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
inline void
SparseMatrixAIJ::addUnsymmetricValuesToSymmetric(const Vector<Int> & is,
                                                 const Vector<Int> & js,
                                                 const Matrix<Real> & values) {
  for (UInt i = 0; i < values.rows(); ++i) {
    UInt c_irn = is(i);
    if (c_irn < size_) {
      for (UInt j = 0; j < values.cols(); ++j) {
        UInt c_jcn = js(j);
        if (c_jcn < size_) {
          if (c_jcn >= c_irn) {
            operator()(c_irn, c_jcn) += values(i, j);
          }
        }
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
inline void
SparseMatrixAIJ::addValuesToUnsymmetric(const Vector<Int> & is,
                                        const Vector<Int> & js,
                                        const Matrix<Real> & values) {
  for (UInt i = 0; i < values.rows(); ++i) {
    UInt c_irn = is(i);
    if (c_irn < size_) {
      for (UInt j = 0; j < values.cols(); ++j) {
        UInt c_jcn = js(j);
        if (c_jcn < size_) {
          operator()(c_irn, c_jcn) += values(i, j);
        }
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
inline void SparseMatrixAIJ::addValues(const Vector<Int> & is,
                                       const Vector<Int> & js,
                                       const Matrix<Real> & values,
                                       MatrixType values_type) {
  if (getMatrixType() == _symmetric)
    if (values_type == _symmetric)
      this->addSymmetricValuesToSymmetric(is, js, values);
    else
      this->addUnsymmetricValuesToSymmetric(is, js, values);
  else
    this->addValuesToUnsymmetric(is, js, values);
}

} // namespace akantu

#endif /* __AKANTU_SPARSE_MATRIX_AIJ_INLINE_IMPL_CC__ */
