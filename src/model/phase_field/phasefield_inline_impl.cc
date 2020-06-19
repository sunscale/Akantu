/* -------------------------------------------------------------------------- */
#include "phase_field_model.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_PHASEFIELD_INLINE_IMPL_CC__
#define __AKANTU_PHASEFIELD_INLINE_IMPL_CC__

namespace akantu {

/* -------------------------------------------------------------------------- */
inline UInt PhaseField::addElement(const ElementType & type, UInt element,
				   const GhostType & ghost_type) {
  Array<UInt> & el_filter = this->element_filter(type, ghost_type);
  el_filter.push_back(element);
  return el_filter.size() - 1;
}

/* -------------------------------------------------------------------------- */
inline UInt PhaseField::addElement(const Element & element) {
  return this->addElement(element.type, element.element, element.ghost_type);
}

  
/* -------------------------------------------------------------------------- */
template <>
inline void PhaseField::registerInternal<Real>(InternalPhaseField<Real> & vect) {
  internal_vectors_real[vect.getID()] = &vect;
}

template <>
inline void PhaseField::registerInternal<UInt>(InternalPhaseField<UInt> & vect) {
  internal_vectors_uint[vect.getID()] = &vect;
}

template <>
inline void PhaseField::registerInternal<bool>(InternalPhaseField<bool> & vect) {
  internal_vectors_bool[vect.getID()] = &vect;
}

/* -------------------------------------------------------------------------- */
template <>
inline void PhaseField::unregisterInternal<Real>(InternalPhaseField<Real> & vect) {
  internal_vectors_real.erase(vect.getID());
}

template <>
inline void PhaseField::unregisterInternal<UInt>(InternalPhaseField<UInt> & vect) {
  internal_vectors_uint.erase(vect.getID());
}

template <>
inline void PhaseField::unregisterInternal<bool>(InternalPhaseField<bool> & vect) {
  internal_vectors_bool.erase(vect.getID());
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline bool PhaseField::isInternal(__attribute__((unused)) const ID & id,
                                 __attribute__((unused))
                                 const ElementKind & element_kind) const {
  AKANTU_TO_IMPLEMENT();
}

template <>
inline bool PhaseField::isInternal<Real>(const ID & id,
                                       const ElementKind & element_kind) const {
  auto internal_array = internal_vectors_real.find(this->getID() + ":" + id);

  if (internal_array == internal_vectors_real.end() ||
      internal_array->second->getElementKind() != element_kind)
    return false;
  return true;
}

  

}

#endif 
