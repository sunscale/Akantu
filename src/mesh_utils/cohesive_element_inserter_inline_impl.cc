/**
 * @file   cohesive_element_inserter_inline_impl.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Wed Nov 05 2014
 * @date last modification: Fri Dec 08 2017
 *
 * @brief  Cohesive element inserter inline functions
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
#include "cohesive_element_inserter.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_COHESIVE_ELEMENT_INSERTER_INLINE_IMPL_CC__
#define __AKANTU_COHESIVE_ELEMENT_INSERTER_INLINE_IMPL_CC__

namespace akantu {

/* -------------------------------------------------------------------------- */
inline UInt
CohesiveElementInserter::getNbData(const Array<Element> & elements,
                                   const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;

  if (tag == _gst_ce_groups) {
    size = elements.size() * (sizeof(bool) + sizeof(unsigned int));
  }

  AKANTU_DEBUG_OUT();
  return size;
}

/* -------------------------------------------------------------------------- */
inline void
CohesiveElementInserter::packData(CommunicationBuffer & buffer,
                                  const Array<Element> & elements,
                                  const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();
  if (tag == _gst_ce_groups)
    packUnpackGroupedInsertionData<true>(buffer, elements);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline void
CohesiveElementInserter::unpackData(CommunicationBuffer & buffer,
                                    const Array<Element> & elements,
                                    const SynchronizationTag & tag) {
  AKANTU_DEBUG_IN();
  if (tag == _gst_ce_groups)
    packUnpackGroupedInsertionData<false>(buffer, elements);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <bool pack_mode>
inline void CohesiveElementInserter::packUnpackGroupedInsertionData(
    CommunicationBuffer & buffer, const Array<Element> & elements) const {

  AKANTU_DEBUG_IN();

  auto current_element_type = _not_defined;
  auto current_ghost_type = _casper;
  auto & physical_names = mesh_facets.registerData<UInt>("physical_names");

  Array<bool> * vect = nullptr;
  Array<UInt> * vect2 = nullptr;

  for (const auto & el : elements) {
    if (el.type != current_element_type ||
        el.ghost_type != current_ghost_type) {
      current_element_type = el.type;
      current_ghost_type = el.ghost_type;
      vect =
          &const_cast<Array<bool> &>(insertion_facets(el.type, el.ghost_type));
      vect2 = &(physical_names(el.type, el.ghost_type));
    }

    Vector<bool> data(vect->storage() + el.element, 1);
    Vector<unsigned int> data2(vect2->storage() + el.element, 1);

    if (pack_mode) {
      buffer << data;
      buffer << data2;
    } else {
      buffer >> data;
      buffer >> data2;
    }
  }

  AKANTU_DEBUG_OUT();
}

} // namespace akantu

#endif /* __AKANTU_COHESIVE_ELEMENT_INSERTER_INLINE_IMPL_CC__ */
