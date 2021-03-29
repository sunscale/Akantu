#include <pybind11/pybind11.h>

#ifndef AKANTU_PY_FRAGMENT_MANAGER_HH_
#define AKANTU_PY_FRAGMENT_MANAGER_HH_

namespace akantu {

void register_fragment_manager(pybind11::module & mod);

} // namespace akantu

#endif // AKANTU_PY_FRAGMENT_MANAGER_HH_
