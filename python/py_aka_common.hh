#include <pybind11/pybind11.h>

#ifndef __AKANTU_PY_AKA_COMMON_HH__
#define __AKANTU_PY_AKA_COMMON_HH__

namespace akantu {

void register_enums(pybind11::module & mod);
void register_initialize(pybind11::module & mod);
void register_functions(pybind11::module & mod);

} // namespace akantu

#endif
