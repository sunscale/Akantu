/**
 * @file   dof_manager_default_inline_impl.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Wed Aug 12 10:57:47 2015
 *
 * @brief  Implementation of the DOFManagerDefault inline functions
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

#ifndef __AKANTU_DOF_MANAGER_DEFAULT_INLINE_IMPL_CC__
#define __AKANTU_DOF_MANAGER_DEFAULT_INLINE_IMPL_CC__

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::addSymmetricElementalMatrixToSymmetric(SparseMatrixAIJ & matrix,
							       const Matrix<Real> & element_mat,
							       const Vector<Real> & equation_numbers,
							       UInt max_size) {
  for (UInt i = 0; i < element_mat.rows(); ++i) {
    UInt c_irn = equation_numbers(i);
    if(c_irn < max_size) {
      for (UInt j = i; j < element_mat.cols(); ++j) {
	UInt c_jcn = equation_numbers(j);
	if(c_jcn < max_size) {
	  matrix(c_irn, c_jcn) += elementary_mat(i, j);
	}
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::addUnsymmetricElementalMatrixToSymmetric(SparseMatrixAIJ & matrix,
								 const Matrix<Real> & element_mat,
								 const Vector<Real> & equation_numbers,
								 UInt max_size) {
  for (UInt i = 0; i < size_mat; ++i) {
    UInt c_irn = equation_numbers(i);
    if(c_irn < size) {
      for (UInt j = 0; j < size_mat; ++j) {
	UInt c_jcn = equation_numbers(j);
	if(c_jcn < size) {
	  if (c_jcn >= c_irn) {
	    matrix(c_irn, c_jcn) += elementary_mat(i, j);
	  }
	}
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::addElementalMatrixToUnsymmetric(SparseMatrixAIJ & matrix,
							const Matrix<Real> & element_mat,
							const Vector<Real> & equation_numbers,
							UInt max_size) {
  for (UInt i = 0; i < element_mat.rows(); ++i) {
    UInt c_irn = equation_numbers(i);
    if(c_irn < max_size) {
      for (UInt j = 0; j < element_mat.cols(); ++j) {
	UInt c_jcn = equation_numbers(j);
	if(c_jcn < max_size) {
	  matrix(c_irn, c_jcn) += elementary_mat(i, j);
	}
      }
    }
  }
}

/* -------------------------------------------------------------------------- */

__END_AKANTU__

#endif /* __AKANTU_DOF_MANAGER_DEFAULT_INLINE_IMPL_CC__ */
