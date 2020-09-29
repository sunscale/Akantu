#include <pybind11/pybind11.h>

#ifndef AKANTU_PY_MATERIAL_HH_
#define AKANTU_PY_MATERIAL_HH_

namespace akantu {

void register_material(pybind11::module & mod);

} // namespace akantu

#endif // AKANTU_PY_MATERIAL_HH_
