#include <pybind11/pybind11.h>

#ifndef AKANTU_PY_HEAT_TRANSFERT_MODEL_HH_
#define AKANTU_PY_HEAT_TRANSFERT_MODEL_HH_

namespace akantu {

void register_heat_transfer_model(pybind11::module & mod);

} // namespace akantu

#endif // AKANTU_PY_HEAT_TRANSFERT_MODEL_HH_
