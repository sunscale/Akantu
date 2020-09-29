#include <pybind11/pybind11.h>

#ifndef __AKANTU_PY_GROUP_MANAGER_HH__
#define __AKANTU_PY_GROUP_MANAGER_HH__

namespace akantu {
void register_group_manager(pybind11::module & mod);
} // namespace akantu

#endif /* __AKANTU_PY_GROUP_MANAGER_HH__ */
