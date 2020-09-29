#include <pybind11/pybind11.h>

#ifndef __AKANTU_PY_MESH_HH__
#define __AKANTU_PY_MESH_HH__

namespace akantu {
void register_mesh(pybind11::module & mod);
} // namespace akantu

#endif // __AKANTU_PY_MESH_HH__
