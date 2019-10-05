/**
 * @file   geometrical_element_property.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Nov 29 2017
 * @date last modification: Thu Nov 30 2017
 *
 * @brief  A Documented file.
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
#include "element_class.hh"
/* -------------------------------------------------------------------------- */
#include <boost/preprocessor.hpp>
/* -------------------------------------------------------------------------- */

namespace akantu {

#define AKANTU_INSTANTIATE_TYPES(r, data, type)                                \
  constexpr std::array<UInt, ElementClass<type>::getNbFacetTypes()>            \
      GeometricalElementProperty<                                              \
          ElementClassProperty<type>::geometrical_type>::nb_facets;            \
  constexpr std::array<UInt, ElementClass<type>::getNbFacetTypes()>            \
      GeometricalElementProperty<                                              \
          ElementClassProperty<type>::geometrical_type>::nb_nodes_per_facet;   \
  constexpr std::array<                                                        \
      UInt, detail::sizeFacetConnectivity<GeometricalElementProperty<          \
                ElementClassProperty<type>::geometrical_type>>()>              \
      GeometricalElementProperty<ElementClassProperty<                         \
          type>::geometrical_type>::facet_connectivity_vect;                   \
  constexpr std::array<ElementType, ElementClass<type>::getNbFacetTypes()>     \
      ElementClassExtraGeometryProperties<type>::facet_type;

BOOST_PP_SEQ_FOR_EACH(AKANTU_INSTANTIATE_TYPES, _,
                      (_not_defined)AKANTU_ek_regular_ELEMENT_TYPE)

#if defined(AKANTU_COHESIVE_ELEMENT)
BOOST_PP_SEQ_FOR_EACH(AKANTU_INSTANTIATE_TYPES, _,
                      AKANTU_ek_cohesive_ELEMENT_TYPE)
#endif

} // namespace akantu
