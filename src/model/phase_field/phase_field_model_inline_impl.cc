/* -------------------------------------------------------------------------- */
#include "aka_named_argument.hh"
#include "phasefield_selector.hh"
#include "phasefield_selector_tmpl.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_PHASE_FIELD_MODEL_INLINE_IMPL_CC__
#define __AKANTU_PHASE_FIELD_MODEL_INLINE_IMPL_CC__

namespace akantu {

/* -------------------------------------------------------------------------- */
inline decltype(auto) PhaseFieldModel::getPhaseFields() {
  return make_dereference_adaptor(phasefields);
}

/* -------------------------------------------------------------------------- */
inline decltype(auto) PhaseFieldModel::getPhaseFields() const {
  return make_dereference_adaptor(phasefields);
}

/* -------------------------------------------------------------------------- */
inline PhaseField & PhaseFieldModel::getPhaseField(UInt mat_index) {
  AKANTU_DEBUG_ASSERT(mat_index < phasefields.size(),
                      "The model " << id << " has no phasefield no "
                                   << mat_index);
  return *phasefields[mat_index];
}

/* -------------------------------------------------------------------------- */
inline const PhaseField & PhaseFieldModel::getPhaseField(UInt mat_index) const {
  AKANTU_DEBUG_ASSERT(mat_index < phasefields.size(),
                      "The model " << id << " has no phasefield no "
                                   << mat_index);
  return *phasefields[mat_index];
}

/* -------------------------------------------------------------------------- */
inline PhaseField & PhaseFieldModel::getPhaseField(const std::string & name) {
  std::map<std::string, UInt>::const_iterator it =
      phasefields_names_to_id.find(name);
  AKANTU_DEBUG_ASSERT(it != phasefields_names_to_id.end(),
                      "The model " << id << " has no phasefield named " << name);
  return *phasefields[it->second];
}

/* -------------------------------------------------------------------------- */
inline UInt
PhaseFieldModel::getPhaseFieldIndex(const std::string & name) const {
  auto it = phasefields_names_to_id.find(name);
  AKANTU_DEBUG_ASSERT(it != phasefields_names_to_id.end(),
                      "The model " << id << " has no phasefield named " << name);
  return it->second;
}

/* -------------------------------------------------------------------------- */
inline const PhaseField &
PhaseFieldModel::getPhaseField(const std::string & name) const {
  auto it = phasefields_names_to_id.find(name);
  AKANTU_DEBUG_ASSERT(it != phasefields_names_to_id.end(),
                      "The model " << id << " has no phasefield named " << name);
  return *phasefields[it->second];
}

/* -------------------------------------------------------------------------- */
} // namespace akantu

#endif /* __AKANTU_PHASE_FIELD_MODEL_INLINE_IMPL_CC__ */
