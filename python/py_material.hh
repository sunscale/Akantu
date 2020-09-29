#include <pybind11/pybind11.h>

#ifndef __AKANTU_PY_MATERIAL_HH__
#define __AKANTU_PY_MATERIAL_HH__

namespace akantu {

void register_material(pybind11::module & mod);

} // namespace akantu

#endif // __AKANTU_PY_MATERIAL_HH__
