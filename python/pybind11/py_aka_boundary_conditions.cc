/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
/* -------------------------------------------------------------------------- */
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
#include "boundary_condition_functor.hh"
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
namespace _aka = akantu;

namespace {

/* -------------------------------------------------------------------------- */

template <typename daughter = _aka::BC::DirichletFunctor>
class PyDirichletFunctor : public daughter {
public:
  /* Inherit the constructors */
  using daughter::daughter;

  /* Trampoline (need one for each virtual function) */
  void operator()(_aka::UInt node, _aka::Vector<bool> & flags,
                  _aka::Vector<_aka::Real> & primal,
                  const _aka::Vector<_aka::Real> & coord) const override {

    PYBIND11_OVERLOAD_NAME(void, daughter, "__call__", operator(), node, flags,
                           primal, coord);
  }
};
/* -------------------------------------------------------------------------- */

template <typename daughter = _aka::BC::NeumannFunctor>
class PyNeumannFunctor : public daughter {
public:
  /* Inherit the constructors */
  using daughter::daughter;

  /* Trampoline (need one for each virtual function) */
  void operator()(const _aka::IntegrationPoint & quad_point,
                  _aka::Vector<_aka::Real> & dual,
                  const _aka::Vector<_aka::Real> & coord,
                  const _aka::Vector<_aka::Real> & normals) const override {

    PYBIND11_OVERLOAD_PURE_NAME(void, daughter, "__call__", operator(),
                                quad_point, dual, coord, normals);
  }
};

/* -------------------------------------------------------------------------- */

template <typename Functor, typename Constructor>
decltype(auto) declareDirichletFunctor(py::module mod, const char * name,
                                       Constructor && cons) {
  py::class_<Functor, PyDirichletFunctor<Functor>, _aka::BC::DirichletFunctor>(
      mod, name)
      .def(cons);
}
template <typename Functor, typename Constructor>
decltype(auto) declareNeumannFunctor(py::module mod, const char * name,
                                     Constructor && cons) {
  py::class_<Functor, PyNeumannFunctor<Functor>, _aka::BC::NeumannFunctor>(mod,
                                                                           name)
      .def(cons);
}

/* -------------------------------------------------------------------------- */

py::module & register_boundary_conditions(py::module & mod) {

  py::class_<_aka::BC::Functor>(mod, "_aka::BCFunctor");
  py::class_<_aka::BC::DirichletFunctor, PyDirichletFunctor<>,
             _aka::BC::Functor>(mod, "DirichletFunctor")
      .def(py::init())
      .def(py::init<_aka::SpatialDirection>());

  py::class_<_aka::BC::NeumannFunctor, PyNeumannFunctor<>, _aka::BC::Functor>(
      mod, "NeumannFunctor")
      .def(py::init());

  declareDirichletFunctor<_aka::BC::FixedValue>(
      mod, "FixedValue", py::init<_aka::Real, _aka::BC::Axis>());

  declareDirichletFunctor<_aka::BC::IncrementValue>(
      mod, "IncrementValue", py::init<_aka::Real, _aka::BC::Axis>());

  declareDirichletFunctor<_aka::BC::Increment>(
      mod, "Increment", py::init<_aka::Vector<_aka::Real> &>());

  declareNeumannFunctor<_aka::BC::FromHigherDim>(
      mod, "FromHigherDim", py::init<_aka::Matrix<_aka::Real> &>());

  declareNeumannFunctor<_aka::BC::FromSameDim>(
      mod, "FromSameDim", py::init<_aka::Vector<_aka::Real> &>());

  declareNeumannFunctor<_aka::BC::FreeBoundary>(mod, "FreeBoundary",
                                                py::init());

  return mod;
}
} // namespace
