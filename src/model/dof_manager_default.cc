/**
 * @file   dof_manager_default.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Tue Aug 11 16:21:01 2015
 *
 * @brief  Implementation of the default DOFManager
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
#include "dof_manager_default.hh"
#include "sparse_matrix_aij.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
inline void DOFManagerDefault::addSymmetricElementalMatrixToSymmetric(
    SparseMatrixAIJ & matrix, const Matrix<Real> & elementary_mat,
    const Vector<UInt> & equation_numbers, UInt max_size) {
  for (UInt i = 0; i < elementary_mat.rows(); ++i) {
    UInt c_irn = equation_numbers(i);
    if (c_irn < max_size) {
      for (UInt j = i; j < elementary_mat.cols(); ++j) {
        UInt c_jcn = equation_numbers(j);
        if (c_jcn < max_size) {
          matrix(c_irn, c_jcn) += elementary_mat(i, j);
        }
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
inline void DOFManagerDefault::addUnsymmetricElementalMatrixToSymmetric(
    SparseMatrixAIJ & matrix, const Matrix<Real> & elementary_mat,
    const Vector<UInt> & equation_numbers, UInt max_size) {
  for (UInt i = 0; i < elementary_mat.rows(); ++i) {
    UInt c_irn = equation_numbers(i);
    if (c_irn < max_size) {
      for (UInt j = 0; j < elementary_mat.cols(); ++j) {
        UInt c_jcn = equation_numbers(j);
        if (c_jcn < max_size) {
          if (c_jcn >= c_irn) {
            matrix(c_irn, c_jcn) += elementary_mat(i, j);
          }
        }
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
inline void DOFManagerDefault::addElementalMatrixToUnsymmetric(
    SparseMatrixAIJ & matrix, const Matrix<Real> & elementary_mat,
    const Vector<UInt> & equation_numbers, UInt max_size) {
  for (UInt i = 0; i < elementary_mat.rows(); ++i) {
    UInt c_irn = equation_numbers(i);
    if (c_irn < max_size) {
      for (UInt j = 0; j < elementary_mat.cols(); ++j) {
        UInt c_jcn = equation_numbers(j);
        if (c_jcn < max_size) {
          matrix(c_irn, c_jcn) += elementary_mat(i, j);
        }
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
DOFManagerDefault::DOFManagerDefault(const Mesh & mesh, const ID & id,
                                     const MemoryID & memory_id)
    : DOFManager(mesh, id, memory_id) {}

/* -------------------------------------------------------------------------- */
DOFManagerDefault::~DOFManagerDefault() {
  AIJMatrixMap::iterator it = this->aij_matrices.begin();
  AIJMatrixMap::iterator end = this->aij_matrices.end();

  for (; it != end; ++it)
    delete it->second;

  this->aij_matrices.clear();
}

/* -------------------------------------------------------------------------- */
SparseMatrix & DOFManagerDefault::getNewMatrix(const ID & matrix_id,
                                               const MatrixType & matrix_type) {
  std::stringstream sstr;
  sstr << this->id << ":" << matrix_id;
  SparseMatrix * sm =
      new SparseMatrixAIJ(*this, matrix_type, sstr.str(), this->memory_id);
  this->registerSparseMatrix(matrix_id, *sm);

  return *sm;
}

/* -------------------------------------------------------------------------- */
SparseMatrix & DOFManagerDefault::getNewMatrix(const ID & matrix_id,
                                               const ID & matrix_to_copy_id) {
  std::stringstream sstr;
  sstr << this->id << ":" << matrix_id;
  SparseMatrixAIJ & sm_to_copy = this->getMatrix(matrix_to_copy_id);
  SparseMatrix * sm =
      new SparseMatrixAIJ(sm_to_copy, sstr.str(), this->memory_id);
  this->registerSparseMatrix(matrix_id, *sm);

  return *sm;
}

/* -------------------------------------------------------------------------- */
SparseMatrixAIJ & DOFManagerDefault::getMatrix(const ID & matrix_id) {
  AIJMatrixMap::iterator it = this->aij_matrices.find(matrix_id);

  if (it == this->aij_matrices.end())
    AKANTU_EXCEPTION("The matrix " << matrix_id
                                   << " does not exists in the DOFManager "
                                   << this->id);

  return *(it->second);
}

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::getSolution(const ID & dof_id,
                                    Array<Real> & solution_array) {
  AKANTU_DEBUG_IN();

  const Array<UInt> & equation_number = this->getLocalEquationNumbers(dof_id);

  UInt nb_degree_of_freedoms =
      solution_array.getSize() * solution_array.getNbComponent();

  AKANTU_DEBUG_ASSERT(
      equation_number.getSize() == nb_degree_of_freedoms,
      "The array to get the solution does not have a correct size."
          << " (" << solution_array.getID() << ")");

  Real * sol_it = solution_array.storage();
  UInt * equ_it = equation_number.storage();

  Array<Real> solution;

  for (UInt d = 0; d < nb_degree_of_freedoms; ++d, ++sol_it, ++equ_it) {
    (*sol_it) = solution(*equ_it);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void
DOFManagerDefault::assembleToResidual(const ID & dof_id,
                                      const Array<Real> & array_to_assemble,
                                      Real scale_factor) {
  AKANTU_DEBUG_IN();

  const Array<UInt> & equation_number = this->getLocalEquationNumbers(dof_id);

  UInt nb_degree_of_freedoms =
      array_to_assemble.getSize() * array_to_assemble.getNbComponent();

  AKANTU_DEBUG_ASSERT(equation_number.getSize() == nb_degree_of_freedoms,
                      "The array to assemble does not have a correct size."
                          << " (" << array_to_assemble.getID() << ")");

  Real * arr_it = array_to_assemble.storage();
  UInt * equ_it = equation_number.storage();

  for (UInt d = 0; d < nb_degree_of_freedoms; ++d, ++arr_it, ++equ_it) {
    residual(*equ_it) += scale_factor * (*arr_it);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::assembleElementalMatricesToMatrix(
    const ID & matrix_id, const ID & dof_id, const Array<Real> & elementary_mat,
    const ElementType & type, const GhostType & ghost_type,
    const MatrixType & elemental_matrix_type,
    const Array<UInt> & filter_elements) {
  AKANTU_DEBUG_IN();

  const Array<UInt> & equation_number = this->getLocalEquationNumbers(dof_id);
  SparseMatrixAIJ & A = this->getMatrix(matrix_id);

  UInt nb_element;
  if (ghost_type == _not_ghost) {
    nb_element = mesh.getNbElement(type);
  } else {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  UInt * filter_it = NULL;
  if (filter_elements != empty_filter) {
    nb_element = filter_elements.getSize();
    filter_it = filter_elements.storage();
  } else {
    nb_element = mesh.getNbElement(type, ghost_type);
  }

  AKANTU_DEBUG_ASSERT(elementary_mat.getSize() == nb_element,
                      "The vector elementary_mat("
                          << elementary_mat.getID()
                          << ") has not the good size.");

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt nb_degree_of_freedom = elementary_mat.getNbComponent() /
                              (nb_nodes_per_element * nb_nodes_per_element);

  const Array<UInt> connectivity = this->mesh.getConnectivity(type, ghost_type);
  Array<UInt>::const_vector_iterator conn_begin =
      connectivity.begin(nb_nodes_per_element);
  Array<UInt>::const_vector_iterator conn_it = conn_begin;

  UInt size_mat = nb_nodes_per_element * nb_degree_of_freedom;

  Vector<UInt> local_eq_nb(nb_degree_of_freedom * nb_nodes_per_element);
  Array<Real>::const_matrix_iterator el_mat_it =
      elementary_mat.begin(size_mat, size_mat);

  for (UInt e = 0; e < nb_element; ++e, ++el_mat_it) {
    if (filter_it != NULL)
      conn_it = conn_begin + *filter_it;

    this->extractElementEquationNumber(equation_number, *conn_it,
                                       nb_degree_of_freedom, local_eq_nb);

    if (filter_it != NULL)
      ++filter_it;
    else
      ++conn_it;

    if (A.getMatrixType() == _symmetric)
      if (elemental_matrix_type == _symmetric)
        this->addSymmetricElementalMatrixToSymmetric(A, *el_mat_it, local_eq_nb,
                                                     A.getSize());
      else
        this->addUnsymmetricElementalMatrixToSymmetric(
            A, *el_mat_it, local_eq_nb, A.getSize());
    else
      this->addElementalMatrixToUnsymmetric(A, *el_mat_it, local_eq_nb,
                                            A.getSize());
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

__END_AKANTU__
