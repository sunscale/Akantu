#include <pybind11/pybind11.h>

#ifndef AKANTU_PY_DUMPABLE_HH_
#define AKANTU_PY_DUMPABLE_HH_

namespace akantu {

void register_dumpable(pybind11::module & mod);

} // namespace akantu

#endif /* AKANTU_PY_DUMPABLE_HH_ */
