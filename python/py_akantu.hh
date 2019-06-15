#include "py_aka_array.hh"

#ifndef __PY_AKANTU_HH__
#define __PY_AKANTU_HH__

namespace pybind11 {
struct module;
} // namespace pybind11

namespace akantu {
[[gnu::visibility("default")]] void register_all(pybind11::module & mod);
}

#endif /* __PY_AKANTU_HH__ */
