/**
 * @file   fe_engine_template_cohesive.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Wed Oct 31 2012
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Specialization of the FEEngineTemplate for cohesive element
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
#include "fe_engine_template.hh"
#include "integrator_gauss.hh"
#include "shape_cohesive.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
/* compatibility functions */
/* -------------------------------------------------------------------------- */
template <>
Real FEEngineTemplate<IntegratorGauss, ShapeLagrange, _ek_cohesive,
                      DefaultIntegrationOrderFunctor>::
    integrate(const Array<Real> & f, const ElementType & type,
              const GhostType & ghost_type,
              const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

#ifndef AKANTU_NDEBUG
  UInt nb_element = mesh.getNbElement(type, ghost_type);
  if (filter_elements != empty_filter)
    nb_element = filter_elements.size();

  UInt nb_quadrature_points = getNbIntegrationPoints(type);

  AKANTU_DEBUG_ASSERT(f.size() == nb_element * nb_quadrature_points,
                      "The vector f(" << f.getID()
                                      << ") has not the good size.");
  AKANTU_DEBUG_ASSERT(f.getNbComponent() == 1,
                      "The vector f("
                          << f.getID()
                          << ") has not the good number of component.");
#endif

  Real integral = 0.;

#define INTEGRATE(type)                                                        \
  integral = integrator.integrate<type>(f, ghost_type, filter_elements);

  AKANTU_BOOST_COHESIVE_ELEMENT_SWITCH(INTEGRATE);
#undef INTEGRATE

  AKANTU_DEBUG_OUT();
  return integral;
}

/* -------------------------------------------------------------------------- */
template <>
void FEEngineTemplate<IntegratorGauss, ShapeLagrange, _ek_cohesive,
                      DefaultIntegrationOrderFunctor>::
    integrate(const Array<Real> & f, Array<Real> & intf,
              UInt nb_degree_of_freedom, const ElementType & type,
              const GhostType & ghost_type,
              const Array<UInt> & filter_elements) const {

#ifndef AKANTU_NDEBUG
  UInt nb_element = mesh.getNbElement(type, ghost_type);
  if (filter_elements != empty_filter)
    nb_element = filter_elements.size();

  UInt nb_quadrature_points = getNbIntegrationPoints(type);

  AKANTU_DEBUG_ASSERT(f.size() == nb_element * nb_quadrature_points,
                      "The vector f(" << f.getID() << " size " << f.size()
                                      << ") has not the good size ("
                                      << nb_element << ").");
  AKANTU_DEBUG_ASSERT(f.getNbComponent() == nb_degree_of_freedom,
                      "The vector f("
                          << f.getID()
                          << ") has not the good number of component.");
  AKANTU_DEBUG_ASSERT(intf.getNbComponent() == nb_degree_of_freedom,
                      "The vector intf("
                          << intf.getID()
                          << ") has not the good number of component.");
  AKANTU_DEBUG_ASSERT(intf.size() == nb_element,
                      "The vector intf(" << intf.getID()
                                         << ") has not the good size.");
#endif

#define INTEGRATE(type)                                                        \
  integrator.integrate<type>(f, intf, nb_degree_of_freedom, ghost_type,        \
                             filter_elements);

  AKANTU_BOOST_COHESIVE_ELEMENT_SWITCH(INTEGRATE);
#undef INTEGRATE
}

/* -------------------------------------------------------------------------- */
template <>
void FEEngineTemplate<IntegratorGauss, ShapeLagrange, _ek_cohesive,
                      DefaultIntegrationOrderFunctor>::
    gradientOnIntegrationPoints(const Array<Real> &, Array<Real> &, const UInt,
                                const ElementType &, const GhostType &,
                                const Array<UInt> &) const {
  AKANTU_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */

} // namespace akantu
