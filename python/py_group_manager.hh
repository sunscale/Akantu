#ifndef __AKANTU_PY_GROUP_MANAGER_HH__
#define __AKANTU_PY_GROUP_MANAGER_HH__

namespace pybind11 {
struct module;
} // namespace pybind11

namespace akantu {
void register_group_manager(pybind11::module & mod);
} // namespace akantu

#endif /* __AKANTU_PY_GROUP_MANAGER_HH__ */
