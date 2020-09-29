/* -------------------------------------------------------------------------- */
#include "py_aka_array.hh"
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
/* -------------------------------------------------------------------------- */

#ifndef PY_AKANTU_HH_
#define PY_AKANTU_HH_

namespace akantu {
void register_all(pybind11::module & mod);
}

#endif /* PY_AKANTU_HH_ */
