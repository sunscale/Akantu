/**
 * @file   fe_engine_template_cohesive.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Wed Oct 31 2012
 * @date last modification: Thu Oct 15 2015
 *
 * @brief  Specialization of the FEEngineTemplate for cohesive element
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include "fe_engine_template.hh"
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
    nb_element = filter_elements.getSize();

  UInt nb_quadrature_points = getNbIntegrationPoints(type);

  AKANTU_DEBUG_ASSERT(f.getSize() == nb_element * nb_quadrature_points,
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
  if (filter_elements == filter_elements)
    nb_element = filter_elements.getSize();

  UInt nb_quadrature_points = getNbIntegrationPoints(type);

  AKANTU_DEBUG_ASSERT(f.getSize() == nb_element * nb_quadrature_points,
                      "The vector f(" << f.getID() << " size " << f.getSize()
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
  AKANTU_DEBUG_ASSERT(intf.getSize() == nb_element,
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
  AKANTU_DEBUG_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
template <>
template <>
void FEEngineTemplate<IntegratorGauss, ShapeLagrange, _ek_cohesive,
                      DefaultIntegrationOrderFunctor>::
    computeNormalsOnIntegrationPoints<_cohesive_1d_2>(
        const Array<Real> &, Array<Real> & normal,
        const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(
      mesh.getSpatialDimension() == 1,
      "Mesh dimension must be 1 to compute normals on 1D cohesive elements!");
  const ElementType type = _cohesive_1d_2;
  const ElementType facet_type = Mesh::getFacetType(type);

  UInt nb_element = mesh.getConnectivity(type, ghost_type).getSize();
  normal.resize(nb_element);

  Array<Element> & facets =
      mesh.getMeshFacets().getSubelementToElement(type, ghost_type);
  Array<std::vector<Element> > & segments =
      mesh.getMeshFacets().getElementToSubelement(facet_type, ghost_type);

  Real values[2];

  for (UInt elem = 0; elem < nb_element; ++elem) {

    for (UInt p = 0; p < 2; ++p) {
      Element f = facets(elem, p);
      Element seg = segments(f.element)[0];
      mesh.getBarycenter(seg.element, seg.type, &(values[p]), seg.ghost_type);
    }

    Real difference = values[0] - values[1];

    AKANTU_DEBUG_ASSERT(difference != 0.,
                        "Error in normal computation for cohesive elements");

    normal(elem) = difference / std::abs(difference);
  }

  AKANTU_DEBUG_OUT();
}

} // akantu
