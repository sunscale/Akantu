#ifndef __AKANTU_PY_AKA_MODEL_HH__
#define __AKANTU_PY_AKA_MODEL_HH__

namespace pybind11 {
struct module;
}

namespace akantu {

void register_model(pybind11::module & mod);

}

#endif
