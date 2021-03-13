#include <pybind11/pybind11.h>

#ifndef AKANTU_PY_AKA_SOLVER_HH_
#define AKANTU_PY_AKA_SOLVER_HH_

namespace akantu {

void register_solvers(pybind11::module & mod);

}

#endif
