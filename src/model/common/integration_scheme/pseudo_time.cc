/**
 * @file   pseudo_time.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Feb 19 2016
 * @date last modification: Wed Jan 31 2018
 *
 * @brief  Implementation of a really simple integration scheme
 *
 * @section LICENSE
 *
 * Copyright (©) 2016-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "pseudo_time.hh"
#include "dof_manager.hh"
#include "sparse_matrix.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
PseudoTime::PseudoTime(DOFManager & dof_manager, const ID & dof_id)
    : IntegrationScheme(dof_manager, dof_id, 0), k_release(0) {}

/* -------------------------------------------------------------------------- */
std::vector<std::string> PseudoTime::getNeededMatrixList() { return {"K"}; }

/* -------------------------------------------------------------------------- */
void PseudoTime::predictor(Real) {}

/* -------------------------------------------------------------------------- */
void PseudoTime::corrector(const SolutionType &, Real) {
  auto & us = this->dof_manager.getDOFs(this->dof_id);
  const auto & deltas = this->dof_manager.getSolution(this->dof_id);
  const auto & blocked_dofs = this->dof_manager.getBlockedDOFs(this->dof_id);

  for (auto && tuple : zip(make_view(us), deltas, make_view(blocked_dofs))) {
    auto & u = std::get<0>(tuple);
    const auto & delta = std::get<1>(tuple);
    const auto & bld = std::get<2>(tuple);
    if (not bld)
      u += delta;
  }
}

/* -------------------------------------------------------------------------- */
void PseudoTime::assembleJacobian(const SolutionType &, Real) {
  SparseMatrix & J = this->dof_manager.getMatrix("J");
  const SparseMatrix & K = this->dof_manager.getMatrix("K");

  if (K.getRelease() == k_release)
    return;

  J.copyProfile(K);
  // J.clear();
  J.add(K);

  k_release = K.getRelease();
}

/* -------------------------------------------------------------------------- */
void PseudoTime::assembleResidual(bool) {}
/* -------------------------------------------------------------------------- */

} // namespace akantu
