#ifndef __AKANTU_PY_AKA_ERROR_HH__
#define __AKANTU_PY_AKA_ERROR_HH__

namespace pybind11 {
struct module;
}

namespace akantu {

void register_error(pybind11::module & mod);

}

#endif
