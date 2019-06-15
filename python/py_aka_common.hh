#ifndef __AKANTU_PY_AKA_COMMON_HH__
#define __AKANTU_PY_AKA_COMMON_HH__

namespace pybind11 {
struct module;
} // namespace pybind11

namespace akantu {

void register_enums(pybind11::module & mod);
void register_initialize(pybind11::module & mod);
void register_functions(pybind11::module & mod);

} // namespace akantu

#endif
