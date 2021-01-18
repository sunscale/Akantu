#include <pybind11/pybind11.h>

#ifndef AKANTU_PY_AKA_ERROR_HH_
#define AKANTU_PY_AKA_ERROR_HH_

namespace akantu {

void register_error(pybind11::module & mod);

}

#endif
