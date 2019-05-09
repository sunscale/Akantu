#ifndef __AKANTU_PY_AKA_MATERIAL_HH__
#define __AKANTU_PY_AKA_MATERIAL_HH__

namespace pybind11 {
struct module;
}

namespace akantu {

void register_material(pybind11::module & mod);

}

#endif
