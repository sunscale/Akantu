#include <pybind11/pybind11.h>

#ifndef AKANTU_PY_SOLID_MECHANICS_MODEL_HH_
#define AKANTU_PY_SOLID_MECHANICS_MODEL_HH_


namespace akantu {

void register_solid_mechanics_model(pybind11::module & mod);

} // namespace akantu

#endif // AKANTU_PY_SOLID_MECHANICS_MODEL_HH_
