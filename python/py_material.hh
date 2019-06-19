#ifndef __AKANTU_PY_MATERIAL_HH__
#define __AKANTU_PY_MATERIAL_HH__

namespace pybind11 {
struct module;
} // namespace pybind11

namespace akantu {

void register_material(pybind11::module & mod);

} // namespace akantu

#endif // __AKANTU_PY_MATERIAL_HH__
