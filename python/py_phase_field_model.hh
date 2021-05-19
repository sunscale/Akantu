#include <pybind11/pybind11.h>

#ifndef __AKANTU_PY_PHASE_FIELD_MODEL_HH__
#define __AKANTU_PY_PHASE_FIELD_MODEL_HH__

namespace akantu {
  void register_phase_field_model(pybind11::module & mod);
  void register_phase_field_coupler(pybind11::module & mod);
} // namespace akantu

#endif // __AKANTU_PY_PHASE_FIELD_MODEL_HH__
