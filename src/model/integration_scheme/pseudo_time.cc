/**
 * @file   pseudo_time.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Wed Feb 17 09:49:10 2016
 *
 * @brief  Implementation of a really simple integration scheme
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
#include "pseudo_time.hh"
#include "dof_manager.hh"
#include "sparse_matrix.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
PseudoTime::PseudoTime(DOFManager & dof_manager, const ID & dof_id)
    : IntegrationScheme(dof_manager, dof_id, 0) {}

/* -------------------------------------------------------------------------- */
void PseudoTime::predictor(__attribute__((unused)) Real delta_t) {}

/* -------------------------------------------------------------------------- */
void PseudoTime::corrector(__attribute__((unused)) const SolutionType & type,
                           __attribute__((unused)) Real delta_t) {
  Array<Real> & u = this->dof_manager.getDOFs(this->dof_id);
  const Array<Real> & delta = this->dof_manager.getSolution(this->dof_id);
  const Array<bool> & blocked_dofs =
      this->dof_manager.getBlockedDOFs(this->dof_id);

  UInt nb_nodes = u.getSize();
  UInt nb_degree_of_freedom = u.getNbComponent() * nb_nodes;

  Array<Real>::scalar_iterator u_it = u.begin_reinterpret(nb_degree_of_freedom);
  Array<Real>::scalar_iterator u_end = u.end_reinterpret(nb_degree_of_freedom);
  Array<Real>::const_scalar_iterator delta_it =
      delta.begin_reinterpret(nb_degree_of_freedom);
  Array<bool>::const_scalar_iterator blocked_dofs_it =
      blocked_dofs.begin_reinterpret(nb_degree_of_freedom);

  for (; u_it != u_end; ++u_it, ++delta_it, ++blocked_dofs_it) {
    if (!(*blocked_dofs_it)) {
      *u_it += *delta_it;
    }
  }
}

/* -------------------------------------------------------------------------- */
void PseudoTime::assembleJacobian(__attribute__((unused)) const SolutionType & type,
                                  __attribute__((unused)) Real delta_t) {
  SparseMatrix & J = this->dof_manager.getMatrix("J");
  const SparseMatrix & K = this->dof_manager.getMatrix("K");

  J.add(K);
}

/* -------------------------------------------------------------------------- */
void PseudoTime::assembleResidual(__attribute__((unused)) bool is_lumped) {}

/* -------------------------------------------------------------------------- */

__END_AKANTU__
