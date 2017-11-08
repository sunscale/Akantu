/**
 * @file   cohesive_element_inserter_inline_impl.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 *
 * @brief  Cohesive element inserter inline functions
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */
/* -------------------------------------------------------------------------- */
#include "cohesive_element_inserter.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_COHESIVE_ELEMENT_INSERTER_INLINE_IMPL_CC__
#define __AKANTU_COHESIVE_ELEMENT_INSERTER_INLINE_IMPL_CC__

namespace akantu {

/* -------------------------------------------------------------------------- */
inline UInt CohesiveElementInserter::getNbData(const Array<Element> & elements,
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
inline void CohesiveElementInserter::packData(CommunicationBuffer & buffer,
                                              const Array<Element> & elements,
                                              const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();
  if (tag == _gst_ce_groups)
    packUnpackGroupedInsertionData<true>(buffer, elements);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline void CohesiveElementInserter::unpackData(CommunicationBuffer & buffer,
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

  ElementType current_element_type = _not_defined;
  GhostType current_ghost_type = _casper;
  ElementTypeMapArray<UInt> & physical_names =
      mesh_facets.registerData<UInt>("physical_names");

  Array<bool> * vect = nullptr;
  Array<unsigned int> * vect2 = nullptr;

  Array<Element>::const_iterator<Element> it = elements.begin();
  Array<Element>::const_iterator<Element> end = elements.end();
  for (; it != end; ++it) {
    const Element & el = *it;
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
