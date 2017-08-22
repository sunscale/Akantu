/**
 * @file   integration_scheme_1st_order.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Fri Oct 23 12:31:32 2015
 *
 * @brief  Implementation of the common functions for 1st order time
 *integrations
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
#include "integration_scheme_1st_order.hh"
#include "dof_manager.hh"
#include "sparse_matrix.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
void IntegrationScheme1stOrder::predictor(Real delta_t) {
  AKANTU_DEBUG_IN();

  Array<Real> & u = this->dof_manager.getDOFs(this->dof_id);
  Array<Real> & u_dot = this->dof_manager.getDOFsDerivatives(this->dof_id, 1);
  const Array<bool> & blocked_dofs =
      this->dof_manager.getBlockedDOFs(this->dof_id);

  this->predictor(delta_t, u, u_dot, blocked_dofs);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void IntegrationScheme1stOrder::corrector(const SolutionType & type,
                                          Real delta_t) {
  AKANTU_DEBUG_IN();

  Array<Real> & u = this->dof_manager.getDOFs(this->dof_id);
  Array<Real> & u_dot = this->dof_manager.getDOFsDerivatives(this->dof_id, 1);

  const Array<Real> & solution = this->dof_manager.getSolution(this->dof_id);
  const Array<bool> & blocked_dofs =
      this->dof_manager.getBlockedDOFs(this->dof_id);

  this->corrector(type, delta_t, u, u_dot, blocked_dofs, solution);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void IntegrationScheme1stOrder::assembleResidual(bool is_lumped) {
  AKANTU_DEBUG_IN();

  const Array<Real> & first_derivative =
      dof_manager.getDOFsDerivatives(this->dof_id, 1);
  if (!is_lumped) {
    this->dof_manager.assembleMatMulVectToResidual(this->dof_id, "M",
                                                   first_derivative, -1);
  } else {
    this->dof_manager.assembleLumpedMatMulVectToResidual(this->dof_id, "M",
                                                         first_derivative, -1);
  }

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */

} // akantu
