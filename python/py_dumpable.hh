#ifndef __AKANTU_PY_DUMPABLE_HH__
#define __AKANTU_PY_DUMPABLE_HH__

namespace pybind11 {
struct module;
} // namespace pybind11

namespace akantu {

void register_dumpable(pybind11::module & mod);

} // namespace akantu

#endif /* __AKANTU_PY_DUMPABLE_HH__ */