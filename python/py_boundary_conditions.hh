#include <pybind11/pybind11.h>

#ifndef __AKANTU_PY_BOUNDARY_CONDITIONS_HH__
#define __AKANTU_PY_BOUNDARY_CONDITIONS_HH__

namespace akantu {

void register_boundary_conditions(pybind11::module & mod);

} // namespace akantu

#endif // __AKANTU_PY_BOUNDARY_CONDITIONS_HH__
