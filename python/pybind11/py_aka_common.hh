#ifndef __AKANTU_PY_AKA_COMMON_HH__
#define __AKANTU_PY_AKA_COMMON_HH__

namespace pybind11 {
struct module;
}

namespace akantu {

void register_enums(pybind11::module & mod);

}

#endif
