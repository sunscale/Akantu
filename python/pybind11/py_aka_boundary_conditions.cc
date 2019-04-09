/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
/* -------------------------------------------------------------------------- */
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
#include "boundary_condition_python_functor.hh"
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
namespace _aka = akantu;

namespace {

/* -------------------------------------------------------------------------- */

class PyDirichletFunctor : public _aka::BC::DirichletFunctor {
public:
  /* Inherit the constructors */
  using _aka::BC::DirichletFunctor::DirichletFunctor;

  /* Trampoline (need one for each virtual function) */
  void operator()(_aka::UInt node, _aka::Vector<bool> & flags,
                  _aka::Vector<_aka::Real> & primal,
                  const _aka::Vector<_aka::Real> & coord) const override {

    PYBIND11_OVERLOAD_NAME(void, _aka::BC::DirichletFunctor,
                           "__call__", operator(), node, flags, primal, coord);
  }
};
/* -------------------------------------------------------------------------- */

py::module & register_boundary_conditions(py::module & mod) {

  py::class_<_aka::BC::Functor>(mod, "BCFunctor");
  py::class_<_aka::BC::DirichletFunctor, PyDirichletFunctor, _aka::BC::Functor>(
      mod, "DirichletFunctor")
      .def(py::init())
      .def(py::init<_aka::SpatialDirection>());


  return mod;
} // namespace

} // namespace
