#include <pybind11/pybind11.h>

#ifndef AKANTU_PY_STRUCTURAL_MECHANICS_MODEL_HH_
#define AKANTU_PY_STRUCTURAL_MECHANICS_MODEL_HH_

namespace akantu {

void register_structural_mechanics_model(pybind11::module & mod);

} // namespace akantu

#endif // AKANTU_PY_STRUCTURAL_MECHANICS_MODEL_HH_
