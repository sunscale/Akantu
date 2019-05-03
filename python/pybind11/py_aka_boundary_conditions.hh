#ifndef __AKANTU_PY_AKA_BOUNDARY_CONDITIONS_HH__
#define __AKANTU_PY_AKA_BOUNDARY_CONDITIONS_HH__

namespace pybind11 {
struct module;
}

namespace akantu {

void register_boundary_conditions(pybind11::module & mod);

}

#endif
