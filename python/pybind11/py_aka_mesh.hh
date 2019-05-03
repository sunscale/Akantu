#ifndef __AKANTU_PY_AKA_MESH_HH__
#define __AKANTU_PY_AKA_MESH_HH__

namespace pybind11 {
struct module;
}

namespace akantu {

void register_mesh(pybind11::module & mod);

}

#endif
