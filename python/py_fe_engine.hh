#include <pybind11/pybind11.h>

#ifndef __AKANTU_PY_FE_ENGINE_HH__
#define __AKANTU_PY_FE_ENGINE_HH__

namespace akantu {

void register_fe_engine(pybind11::module & mod);

} // namespace akantu

#endif // __AKANTU_PY_FE_ENGINE_HH__
