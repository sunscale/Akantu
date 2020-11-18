#include <pybind11/pybind11.h>

#ifndef AKANTU_PY_FE_ENGINE_HH_
#define AKANTU_PY_FE_ENGINE_HH_

namespace akantu {

void register_fe_engine(pybind11::module & mod);

} // namespace akantu

#endif // AKANTU_PY_FE_ENGINE_HH_
