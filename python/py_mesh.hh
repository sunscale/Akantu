#include <pybind11/pybind11.h>

#ifndef AKANTU_PY_MESH_HH_
#define AKANTU_PY_MESH_HH_

namespace akantu {
void register_mesh(pybind11::module & mod);
} // namespace akantu

#endif // AKANTU_PY_MESH_HH_
