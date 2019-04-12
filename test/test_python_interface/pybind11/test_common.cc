#include <aka_array.hh>
#include <pybind11/pybind11.h>

#include "../../python/pybind11/py_aka_array.cc"
#include "../../python/pybind11/py_aka_boundary_conditions.cc"
#include "../../python/pybind11/py_aka_common.cc"
#include "../../python/pybind11/py_aka_solid_mechanics_model.cc"

namespace py = pybind11;
namespace _aka = akantu;

std::map<long, std::shared_ptr<_aka::Array<_aka::Real>>> arrays;

PYBIND11_MODULE(py11_akantu_test_common, mod) {
  mod.doc() = "Akantu Test function for common ";

  register_enums(mod);
  register_boundary_conditions(mod);
  register_solid_mechanics_models(mod);

  mod.def("createData",
          [&](_aka::UInt size, _aka::UInt nb_components) {
            auto ptr =
                std::make_shared<_aka::Array<_aka::Real>>(size, nb_components);
            ptr->clear();
            long addr = (long)ptr->storage();
            py::print("initial pointer: " + std::to_string(addr));
            arrays[addr] = ptr;
            return std::tuple<long, _aka::Array<_aka::Real> &>(addr, *ptr);
          },
          py::return_value_policy::reference);
  mod.def("getData",
          [&](long addr) -> _aka::Array<_aka::Real> & {
            auto & array = *arrays[addr];
            py::print("gotten pointer: " +
                      std::to_string((long)array.storage()));
            return array;
          },
          py::return_value_policy::reference);

  mod.def("getCopyData",
          [&](long addr) -> _aka::Array<_aka::Real> {
            auto & array = *arrays[addr];
            py::print("gotten pointer: " +
                      std::to_string((long)array.storage()));
            return array;
          },
          py::return_value_policy::copy);

  mod.def("getRawPointer", [](_aka::Array<_aka::Real> & _data) {
    py::print("received proxy: " + std::to_string((long)&_data));
    py::print("raw pointer: " + std::to_string((long)_data.storage()));
    return (long)_data.storage();
  });
} // Module akantu_test_common
