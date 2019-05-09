#ifndef __AKANTU_PY_AKA_FE_ENGINE_HH__
#define __AKANTU_PY_AKA_FE_ENGINE_HH__

namespace pybind11 {
struct module;
}

namespace akantu {

void register_fe_engine(pybind11::module & mod);

}

#endif
