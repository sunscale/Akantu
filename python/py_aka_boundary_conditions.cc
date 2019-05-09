/* -------------------------------------------------------------------------- */
#include "py_aka_boundary_conditions.hh"
#include "py_aka_array.hh"
/* -------------------------------------------------------------------------- */
#include <aka_common.hh>
#include <boundary_condition_functor.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;

namespace akantu {

/* -------------------------------------------------------------------------- */

template <typename daughter = BC::Dirichlet::DirichletFunctor>
class PyDirichletFunctor : public daughter {
public:
  /* Inherit the constructors */
  using daughter::daughter;

  /* Trampoline (need one for each virtual function) */
  void operator()(UInt node, Vector<bool> & flags, Vector<Real> & primal,
                  const Vector<Real> & coord) const override {

    PYBIND11_OVERLOAD_NAME(void, daughter, "__call__", operator(), node, flags,
                           primal, coord);
  }
};
/* -------------------------------------------------------------------------- */

template <typename daughter = BC::Neumann::NeumannFunctor>
class PyNeumannFunctor : public daughter {
public:
  /* Inherit the constructors */
  using daughter::daughter;

  /* Trampoline (need one for each virtual function) */
  void operator()(const IntegrationPoint & quad_point, Vector<Real> & dual,
                  const Vector<Real> & coord,
                  const Vector<Real> & normals) const override {

    PYBIND11_OVERLOAD_PURE_NAME(void, daughter, "__call__", operator(),
                                quad_point, dual, coord, normals);
  }
};

/* -------------------------------------------------------------------------- */

template <typename Functor, typename Constructor>
decltype(auto) declareDirichletFunctor(py::module mod, const char * name,
                                       Constructor && cons) {
  py::class_<Functor, PyDirichletFunctor<Functor>,
             BC::Dirichlet::DirichletFunctor>(mod, name)
      .def(cons);
}
template <typename Functor, typename Constructor>
decltype(auto) declareNeumannFunctor(py::module mod, const char * name,
                                     Constructor && cons) {
  py::class_<Functor, PyNeumannFunctor<Functor>, BC::Neumann::NeumannFunctor>(
      mod, name)
      .def(cons);
}

/* -------------------------------------------------------------------------- */
__attribute__((visibility("default"))) void
register_boundary_conditions(py::module & mod) {

  py::class_<BC::Functor>(mod, "BCFunctor");
  py::class_<BC::Dirichlet::DirichletFunctor, PyDirichletFunctor<>,
             BC::Functor>(mod, "DirichletFunctor")
      .def(py::init())
      .def(py::init<SpatialDirection>());

  py::class_<BC::Neumann::NeumannFunctor, PyNeumannFunctor<>, BC::Functor>(
      mod, "NeumannFunctor")
      .def(py::init());

  declareDirichletFunctor<BC::Dirichlet::FixedValue>(
      mod, "FixedValue", py::init<Real, BC::Axis>());

  declareDirichletFunctor<BC::Dirichlet::IncrementValue>(
      mod, "IncrementValue", py::init<Real, BC::Axis>());

  declareDirichletFunctor<BC::Dirichlet::Increment>(mod, "Increment",
                                                    py::init<Vector<Real> &>());

  declareNeumannFunctor<BC::Neumann::FromHigherDim>(mod, "FromHigherDim",
                                                    py::init<Matrix<Real> &>());

  declareNeumannFunctor<BC::Neumann::FromSameDim>(mod, "FromSameDim",
                                                  py::init<Vector<Real> &>());

  declareNeumannFunctor<BC::Neumann::FreeBoundary>(mod, "FreeBoundary",
                                                   py::init());
}
} // namespace akantu
