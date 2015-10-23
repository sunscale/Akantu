/**
 * @file   integration_scheme_2nd_order.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Tue Oct 20 10:41:12 2015
 *
 * @brief  Implementation of the common part of 2nd order integration schemes
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
#include "integration_scheme_2nd_order.hh"
#include "dof_manager.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
void IntegrationScheme2ndOrder::predictor(const ID & dof_id,
                                          Real delta_t) {
  AKANTU_DEBUG_IN();

  Array<Real> & u = this->dof_manager.getDOFs(dof_id);
  Array<Real> & u_dot = this->dof_manager.getDOFsDerivatives(dof_id, 1);
  Array<Real> & u_dot_dot = this->dof_manager.getDOFsDerivatives(dof_id, 2);
  const Array<bool> & blocked_dofs = this->dof_manager.getBlockedDOFs(dof_id);

  this->predictor(delta_t, u, u_dot, u_dot_dot, blocked_dofs);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void IntegrationScheme2ndOrder::corrector(const SolutionType & type,
                                          const ID & dof_id,
                                          Real delta_t) {
  AKANTU_DEBUG_IN();

  Array<Real> & u = this->dof_manager.getDOFs(dof_id);
  Array<Real> & u_dot = this->dof_manager.getDOFsDerivatives(dof_id, 1);
  Array<Real> & u_dot_dot = this->dof_manager.getDOFsDerivatives(dof_id, 2);

  const Array<Real> & solution = this->dof_manager.getSolution(dof_id);
  const Array<bool> & blocked_dofs = this->dof_manager.getBlockedDOFs(dof_id);

  this->corrector(type, delta_t, u, u_dot, u_dot_dot, blocked_dofs, solution);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

__END_AKANTU__
