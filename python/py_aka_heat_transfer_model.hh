#ifndef __AKANTU_PY_AKA_HEAT_TRANSFERT_MODEL_HH__
#define __AKANTU_PY_AKA_HEAT_TRANSFERT_MODEL_HH__

namespace pybind11 {
struct module;
}

namespace akantu {

void register_heat_transfer_model(pybind11::module & mod);

}

#endif
