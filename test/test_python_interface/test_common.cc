/* -------------------------------------------------------------------------- */
#include "py_akantu.hh"
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
/* -------------------------------------------------------------------------- */
#include <map>
/* -------------------------------------------------------------------------- */

namespace py = pybind11;
namespace _aka = akantu;

std::map<long, std::shared_ptr<_aka::Array<_aka::Real>>> arrays;
std::map<long, std::shared_ptr<_aka::Vector<_aka::Real>>> vectors;
std::map<long, std::shared_ptr<_aka::Matrix<_aka::Real>>> matrices;

PYBIND11_MODULE(py11_akantu_test_common, mod) {
  mod.doc() = "Akantu Test function for common ";

  mod.def("createArray",
          [&](_aka::UInt size, _aka::UInt nb_components) {
            auto ptr =
                std::make_shared<_aka::Array<_aka::Real>>(size, nb_components);
            ptr->zero();
            long addr = (long)ptr->storage();
            py::print("initial pointer: " + std::to_string(addr));
            arrays[addr] = ptr;
            return std::tuple<long, _aka::Array<_aka::Real> &>(addr, *ptr);
          },
          py::return_value_policy::reference);
  mod.def("getArray",
          [&](long addr) -> _aka::Array<_aka::Real> & {
            auto & array = *arrays[addr];
            py::print("gotten pointer: " +
                      std::to_string((long)array.storage()));
            return array;
          },
          py::return_value_policy::reference);

  mod.def("copyArray",
          [&](long addr) -> _aka::Array<_aka::Real> {
            auto & array = *arrays[addr];
            py::print("gotten pointer: " +
                      std::to_string((long)array.storage()));
            return array;
          },
          py::return_value_policy::copy);

  mod.def("getRawPointerArray", [](_aka::Array<_aka::Real> & _data) {
    py::print("received proxy: " + std::to_string((long)&_data));
    py::print("raw pointer: " + std::to_string((long)_data.storage()));
    return (long)_data.storage();
  });

  mod.def("createVector",
          [&](_aka::UInt size) {
            auto ptr = std::make_shared<_aka::Vector<_aka::Real>>(size);
            ptr->zero();
            long addr = (long)ptr->storage();
            py::print("initial pointer: " + std::to_string(addr));
            vectors[addr] = ptr;
            return std::tuple<long, _aka::Vector<_aka::Real> &>(addr, *ptr);
          },
          py::return_value_policy::reference);
  mod.def("getVector",
          [&](long addr) -> _aka::Vector<_aka::Real> & {
            auto & vector = *vectors[addr];
            py::print("gotten pointer: " +
                      std::to_string((long)vector.storage()));
            return vector;
          },
          py::return_value_policy::reference);

  mod.def("copyVector",
          [&](long addr) -> _aka::Vector<_aka::Real> {
            auto & vector = *vectors[addr];
            py::print("gotten pointer: " +
                      std::to_string((long)vector.storage()));
            return vector;
          },
          py::return_value_policy::copy);

  mod.def("getRawPointerVector", [](_aka::Vector<_aka::Real> & _data) {
    py::print("received proxy: " + std::to_string((long)&_data));
    py::print("raw pointer: " + std::to_string((long)_data.storage()));
    return (long)_data.storage();
  });

  mod.def("createMatrix",
          [&](_aka::UInt size1, _aka::UInt size2) {
            auto ptr = std::make_shared<_aka::Matrix<_aka::Real>>(size1, size2);
            ptr->zero();
            long addr = (long)ptr->storage();
            py::print("initial pointer: " + std::to_string(addr));
            matrices[addr] = ptr;
            return std::tuple<long, _aka::Matrix<_aka::Real> &>(addr, *ptr);
          },
          py::return_value_policy::reference);
  mod.def("getMatrix",
          [&](long addr) -> _aka::Matrix<_aka::Real> & {
            auto & matrix = *matrices[addr];
            py::print("gotten pointer: " +
                      std::to_string((long)matrix.storage()));
            return matrix;
          },
          py::return_value_policy::reference);

  mod.def("copyMatrix",
          [&](long addr) -> _aka::Matrix<_aka::Real> {
            auto & matrix = *matrices[addr];
            py::print("gotten pointer: " +
                      std::to_string((long)matrix.storage()));
            return matrix;
          },
          py::return_value_policy::copy);

  mod.def("getRawPointerMatrix", [](_aka::Matrix<_aka::Real> & _data) {
    py::print("received proxy: " + std::to_string((long)&_data));
    py::print("raw pointer: " + std::to_string((long)_data.storage()));
    return (long)_data.storage();
  });
} // Module akantu_test_common
