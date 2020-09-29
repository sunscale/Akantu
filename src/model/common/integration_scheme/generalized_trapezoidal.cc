/**
 * @file   generalized_trapezoidal.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Oct 23 2015
 * @date last modification: Wed Jan 31 2018
 *
 * @brief  implementation of inline functions
 *
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
#include "generalized_trapezoidal.hh"
#include "aka_array.hh"
#include "dof_manager.hh"
#include "mesh.hh"
#include "sparse_matrix.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
GeneralizedTrapezoidal::GeneralizedTrapezoidal(DOFManager & dof_manager,
                                               const ID & dof_id, Real alpha)
    : IntegrationScheme1stOrder(dof_manager, dof_id), alpha(alpha) {

  this->registerParam("alpha", this->alpha, alpha, _pat_parsmod,
                      "The alpha parameter");
}

/* -------------------------------------------------------------------------- */
void GeneralizedTrapezoidal::predictor(Real delta_t, Array<Real> & u,
                                       Array<Real> & u_dot,
                                       const Array<bool> & blocked_dofs) const {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = u.size();
  UInt nb_degree_of_freedom = u.getNbComponent() * nb_nodes;

  Real * u_val = u.storage();
  Real * u_dot_val = u_dot.storage();
  bool * blocked_dofs_val = blocked_dofs.storage();

  for (UInt d = 0; d < nb_degree_of_freedom; d++) {
    if (!(*blocked_dofs_val)) {
      *u_val += (1. - alpha) * delta_t * *u_dot_val;
    }
    u_val++;
    u_dot_val++;
    blocked_dofs_val++;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void GeneralizedTrapezoidal::corrector(const SolutionType & type, Real delta_t,
                                       Array<Real> & u, Array<Real> & u_dot,
                                       const Array<bool> & blocked_dofs,
                                       const Array<Real> & delta) const {
  AKANTU_DEBUG_IN();

  switch (type) {
  case _temperature:
    this->allCorrector<_temperature>(delta_t, u, u_dot, blocked_dofs, delta);
    break;
  case _temperature_rate:
    this->allCorrector<_temperature_rate>(delta_t, u, u_dot, blocked_dofs,
                                          delta);
    break;
  default:
    AKANTU_EXCEPTION("The corrector type : "
                     << type
                     << " is not supported by this type of integration scheme");
  }

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
Real GeneralizedTrapezoidal::getTemperatureCoefficient(
    const SolutionType & type, Real delta_t) const {
  switch (type) {
  case _temperature:
    return 1.;
  case _temperature_rate:
    return alpha * delta_t;
  default:
    AKANTU_EXCEPTION("The corrector type : "
                     << type
                     << " is not supported by this type of integration scheme");
  }
}

/* -------------------------------------------------------------------------- */
Real GeneralizedTrapezoidal::getTemperatureRateCoefficient(
    const SolutionType & type, Real delta_t) const {
  switch (type) {
  case _temperature:
    return 1. / (alpha * delta_t);
  case _temperature_rate:
    return 1.;
  default:
    AKANTU_EXCEPTION("The corrector type : "
                     << type
                     << " is not supported by this type of integration scheme");
  }
}

/* -------------------------------------------------------------------------- */
template <IntegrationScheme::SolutionType type>
void GeneralizedTrapezoidal::allCorrector(Real delta_t, Array<Real> & u,
                                          Array<Real> & u_dot,
                                          const Array<bool> & blocked_dofs,
                                          const Array<Real> & delta) const {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = u.size();
  UInt nb_degree_of_freedom = u.getNbComponent() * nb_nodes;

  Real e = getTemperatureCoefficient(type, delta_t);
  Real d = getTemperatureRateCoefficient(type, delta_t);

  Real * u_val = u.storage();
  Real * u_dot_val = u_dot.storage();
  Real * delta_val = delta.storage();
  bool * blocked_dofs_val = blocked_dofs.storage();

  for (UInt dof = 0; dof < nb_degree_of_freedom; dof++) {
    if (!(*blocked_dofs_val)) {
      *u_val += e * *delta_val;
      *u_dot_val += d * *delta_val;
    }
    u_val++;
    u_dot_val++;
    delta_val++;
    blocked_dofs_val++;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void GeneralizedTrapezoidal::assembleJacobian(const SolutionType & type,
                                              Real delta_t) {
  AKANTU_DEBUG_IN();

  SparseMatrix & J = this->dof_manager.getMatrix("J");

  const SparseMatrix & M = this->dof_manager.getMatrix("M");
  const SparseMatrix & K = this->dof_manager.getMatrix("K");

  bool does_j_need_update = false;
  does_j_need_update |= M.getRelease() != m_release;
  does_j_need_update |= K.getRelease() != k_release;
  if (!does_j_need_update) {
    AKANTU_DEBUG_OUT();
    return;
  }

  J.copyProfile(K);
  // J.zero();

  Real c = this->getTemperatureRateCoefficient(type, delta_t);
  Real e = this->getTemperatureCoefficient(type, delta_t);

  J.add(M, e);
  J.add(K, c);

  m_release = M.getRelease();
  k_release = K.getRelease();

  AKANTU_DEBUG_OUT();
}

} // namespace akantu
