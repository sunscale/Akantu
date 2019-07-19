/**
 * @file   sparse_matrix.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Dec 13 2010
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  implementation of the SparseMatrix class
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
#include <fstream>
/* -------------------------------------------------------------------------- */
#include "communicator.hh"
#include "dof_manager.hh"
#include "sparse_matrix.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
SparseMatrix::SparseMatrix(DOFManager & dof_manager,
                           const MatrixType & matrix_type, const ID & id)
    : id(id), _dof_manager(dof_manager), matrix_type(matrix_type),
      size_(dof_manager.getSystemSize()), nb_non_zero(0) {
  AKANTU_DEBUG_IN();

  const auto & comm = _dof_manager.getCommunicator();
  this->nb_proc = comm.getNbProc();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
SparseMatrix::SparseMatrix(const SparseMatrix & matrix, const ID & id)
    : SparseMatrix(matrix._dof_manager, matrix.matrix_type, id) {
  nb_non_zero = matrix.nb_non_zero;
}

/* -------------------------------------------------------------------------- */
SparseMatrix::~SparseMatrix() = default;

// /* --------------------------------------------------------------------------
// */ Array<Real> & operator*=(SolverVector & vect, const SparseMatrix & mat) {
//   Array<Real> tmp(vect.size(), vect.getNbComponent(), 0.);
//   mat.matVecMul(vect, tmp);

//   vect.copy(tmp);
//   return vect;
// }

/* -------------------------------------------------------------------------- */
void SparseMatrix::add(const SparseMatrix & B, Real alpha) {
  B.addMeTo(*this, alpha);
}

/* -------------------------------------------------------------------------- */

} // namespace akantu
