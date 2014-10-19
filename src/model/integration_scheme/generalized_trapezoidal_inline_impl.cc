/**
 * @file   generalized_trapezoidal_inline_impl.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Jul 04 2011
 * @date last modification: Thu Jun 05 2014
 *
 * @brief  implementation of inline functions
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

inline void GeneralizedTrapezoidal::integrationSchemePred(Real delta_t,
							  Array<Real> & u,
							  Array<Real> & u_dot,
							  Array<bool> & blocked_dofs) {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = u.getSize();
  UInt nb_degree_of_freedom = u.getNbComponent() * nb_nodes;

  Real * u_val         = u.storage();
  Real * u_dot_val     = u_dot.storage();
  bool * blocked_dofs_val  = blocked_dofs.storage();

  for (UInt d = 0; d < nb_degree_of_freedom; d++) {
    if(!(*blocked_dofs_val)) {
      *u_val += delta_t * *u_dot_val;
    }
    u_val++;
    u_dot_val++;
    blocked_dofs_val++;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline void GeneralizedTrapezoidal::integrationSchemeCorrTemp(Real delta_t,
						       Array<Real> & u,
						       Array<Real> & u_dot,
						       Array<bool> & blocked_dofs,
						       Array<Real> & delta) {
  AKANTU_DEBUG_IN();

  integrationSchemeCorr<GeneralizedTrapezoidal::_temperature_corrector>(delta_t,
									u,
									u_dot,
									blocked_dofs,
									delta);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline void GeneralizedTrapezoidal::integrationSchemeCorrTempRate(Real delta_t,
							   Array<Real> & u,
							   Array<Real> & u_dot,
							   Array<bool> & blocked_dofs,
							   Array<Real> & delta) {
  AKANTU_DEBUG_IN();

  integrationSchemeCorr<GeneralizedTrapezoidal::_temperature_rate_corrector>(delta_t,
									     u,
									     u_dot,
									     blocked_dofs,
									     delta);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<>
inline Real GeneralizedTrapezoidal::getTemperatureCoefficient<GeneralizedTrapezoidal::_temperature_corrector>(__attribute__ ((unused)) Real delta_t) {
  return 1.;
}

template<>
inline Real GeneralizedTrapezoidal::getTemperatureRateCoefficient<GeneralizedTrapezoidal::_temperature_corrector>(Real delta_t) {
  return 1./(alpha * delta_t);
}

/* -------------------------------------------------------------------------- */
template<>
inline Real GeneralizedTrapezoidal::getTemperatureCoefficient<GeneralizedTrapezoidal::_temperature_rate_corrector>(Real delta_t) {
  return alpha * delta_t;
}

template<>
inline Real GeneralizedTrapezoidal::getTemperatureRateCoefficient<GeneralizedTrapezoidal::_temperature_rate_corrector>(__attribute__ ((unused)) Real delta_t) {
  return 1.;
}

/* -------------------------------------------------------------------------- */
template<GeneralizedTrapezoidal::IntegrationSchemeCorrectorType type>
inline void GeneralizedTrapezoidal::integrationSchemeCorr(Real delta_t,
						   Array<Real> & u,
						   Array<Real> & u_dot,
						   Array<bool> & blocked_dofs,
						   Array<Real> & delta) {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = u.getSize();
  UInt nb_degree_of_freedom = u.getNbComponent() * nb_nodes;

  Real e = getTemperatureCoefficient<type>(delta_t);
  Real d = getTemperatureRateCoefficient<type>(delta_t);

  Real * u_val         = u.storage();
  Real * u_dot_val     = u_dot.storage();
  Real * delta_val     = delta.storage();
  bool * blocked_dofs_val  = blocked_dofs.storage();

  for (UInt dof = 0; dof < nb_degree_of_freedom; dof++) {
    if(!(*blocked_dofs_val)) {
      *u_val         += e * *delta_val;
      *u_dot_val     += d * *delta_val;
    }
    u_val++;
    u_dot_val++;
    delta_val++;
    blocked_dofs_val++;
  }

  AKANTU_DEBUG_OUT();
}
