#include <pybind11/pybind11.h>

#ifndef AKANTU_PY_AKA_COMMON_HH_
#define AKANTU_PY_AKA_COMMON_HH_

namespace akantu {

void register_enums(pybind11::module & mod);
void register_initialize(pybind11::module & mod);
void register_functions(pybind11::module & mod);

} // namespace akantu

#endif
