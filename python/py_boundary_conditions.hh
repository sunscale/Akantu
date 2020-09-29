#include <pybind11/pybind11.h>

#ifndef AKANTU_PY_BOUNDARY_CONDITIONS_HH_
#define AKANTU_PY_BOUNDARY_CONDITIONS_HH_

namespace akantu {

void register_boundary_conditions(pybind11::module & mod);

} // namespace akantu

#endif // AKANTU_PY_BOUNDARY_CONDITIONS_HH_
