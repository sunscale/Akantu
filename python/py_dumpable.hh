#include <pybind11/pybind11.h>

#ifndef __AKANTU_PY_DUMPABLE_HH__
#define __AKANTU_PY_DUMPABLE_HH__

namespace akantu {

void register_dumpable(pybind11::module & mod);

} // namespace akantu

#endif /* __AKANTU_PY_DUMPABLE_HH__ */
