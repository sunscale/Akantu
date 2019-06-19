#ifndef __AKANTU_PY_FE_ENGINE_HH__
#define __AKANTU_PY_FE_ENGINE_HH__

namespace pybind11 {
struct module;
} // namespace pybind11

namespace akantu {

void register_fe_engine(pybind11::module & mod);

} // namespace akantu

#endif // __AKANTU_PY_FE_ENGINE_HH__
