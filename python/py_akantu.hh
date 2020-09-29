/* -------------------------------------------------------------------------- */
#include "py_aka_array.hh"
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
/* -------------------------------------------------------------------------- */

#ifndef __PY_AKANTU_HH__
#define __PY_AKANTU_HH__

namespace akantu {
void register_all(pybind11::module & mod);
}

#endif /* __PY_AKANTU_HH__ */
