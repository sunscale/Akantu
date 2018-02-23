/**
 * @file   integration_scheme_2nd_order.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Oct 23 2015
 * @date last modification: Wed Jan 31 2018
 *
 * @brief  Implementation of the common part of 2nd order integration schemes
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
#include "integration_scheme_2nd_order.hh"
#include "dof_manager.hh"
#include "sparse_matrix.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
std::vector<std::string> IntegrationScheme2ndOrder::getNeededMatrixList() {
  return {"K", "M", "C"};
}

/* -------------------------------------------------------------------------- */
void IntegrationScheme2ndOrder::predictor(Real delta_t) {
  AKANTU_DEBUG_IN();

  Array<Real> & u = this->dof_manager.getDOFs(this->dof_id);
  Array<Real> & u_dot = this->dof_manager.getDOFsDerivatives(this->dof_id, 1);
  Array<Real> & u_dot_dot =
      this->dof_manager.getDOFsDerivatives(this->dof_id, 2);
  const Array<bool> & blocked_dofs =
      this->dof_manager.getBlockedDOFs(this->dof_id);

  this->predictor(delta_t, u, u_dot, u_dot_dot, blocked_dofs);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void IntegrationScheme2ndOrder::corrector(const SolutionType & type,
                                          Real delta_t) {
  AKANTU_DEBUG_IN();

  Array<Real> & u = this->dof_manager.getDOFs(this->dof_id);
  Array<Real> & u_dot = this->dof_manager.getDOFsDerivatives(this->dof_id, 1);
  Array<Real> & u_dot_dot =
      this->dof_manager.getDOFsDerivatives(this->dof_id, 2);

  const Array<Real> & solution = this->dof_manager.getSolution(this->dof_id);
  const Array<bool> & blocked_dofs =
      this->dof_manager.getBlockedDOFs(this->dof_id);

  this->corrector(type, delta_t, u, u_dot, u_dot_dot, blocked_dofs, solution);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void IntegrationScheme2ndOrder::assembleResidual(bool is_lumped) {
  AKANTU_DEBUG_IN();

  if (this->dof_manager.hasMatrix("C")) {
    const Array<Real> & first_derivative =
        this->dof_manager.getDOFsDerivatives(this->dof_id, 1);

    this->dof_manager.assembleMatMulVectToResidual(this->dof_id, "C",
                                                   first_derivative, -1);
  }

  const Array<Real> & second_derivative =
      this->dof_manager.getDOFsDerivatives(this->dof_id, 2);

  if (not is_lumped) {
    this->dof_manager.assembleMatMulVectToResidual(this->dof_id, "M",
                                                   second_derivative, -1);
  } else {
    this->dof_manager.assembleLumpedMatMulVectToResidual(this->dof_id, "M",
                                                         second_derivative, -1);
  }

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */

} // namespace akantu
