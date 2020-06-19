/**
 * @file   material_cohesive_inline_impl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Wed Aug 04 2010
 * @date last modification: Mon Feb 19 2018
 *
 * @brief  MaterialCohesive inline implementation
 *
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
#include "material_cohesive.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
inline UInt MaterialCohesive::addFacet(const Element & element) {
  Array<UInt> & f_filter = facet_filter(element.type, element.ghost_type);
  f_filter.push_back(element.element);
  return f_filter.size() - 1;
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void MaterialCohesive::computeNormal(const Array<Real> & /*position*/,
                                     Array<Real> & /*normal*/,
                                     GhostType /*ghost_type*/) {}

/* -------------------------------------------------------------------------- */
inline UInt MaterialCohesive::getNbData(const Array<Element> & elements,
                                        const SynchronizationTag & tag) const {

  switch (tag) {
  case SynchronizationTag::_smm_stress: {
    return 2 * spatial_dimension * sizeof(Real) *
           this->getModel().getNbIntegrationPoints(elements,
                                                   "CohesiveFEEngine");
  }
  case SynchronizationTag::_smmc_damage: {
    return sizeof(Real) * this->getModel().getNbIntegrationPoints(
                              elements, "CohesiveFEEngine");
  }
  default: {
  }
  }

  return 0;
}

/* -------------------------------------------------------------------------- */
inline void MaterialCohesive::packData(CommunicationBuffer & buffer,
                                       const Array<Element> & elements,
                                       const SynchronizationTag & tag) const {
  switch (tag) {
  case SynchronizationTag::_smm_stress: {
    packElementDataHelper(tractions, buffer, elements, "CohesiveFEEngine");
    packElementDataHelper(contact_tractions, buffer, elements,
                          "CohesiveFEEngine");
    break;
  }
  case SynchronizationTag::_smmc_damage:
    packElementDataHelper(damage, buffer, elements, "CohesiveFEEngine");
    break;
  default: {
  }
  }
}

/* -------------------------------------------------------------------------- */
inline void MaterialCohesive::unpackData(CommunicationBuffer & buffer,
                                         const Array<Element> & elements,
                                         const SynchronizationTag & tag) {
  switch (tag) {
  case SynchronizationTag::_smm_stress: {
    unpackElementDataHelper(tractions, buffer, elements, "CohesiveFEEngine");
    unpackElementDataHelper(contact_tractions, buffer, elements,
                            "CohesiveFEEngine");
    break;
  }
  case SynchronizationTag::_smmc_damage:
    unpackElementDataHelper(damage, buffer, elements, "CohesiveFEEngine");
    break;
  default: {
  }
  }
}
} // namespace akantu
