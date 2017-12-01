/**
 * @file   test_solid_mechanics_model_cube3d.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date creation: Wed Aug 04 2010
 * @date last modification: Thu Aug 06 2015
 *
 * @brief  test of the class SolidMechanicsModel on the 3d cube
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "boundary_condition_functor.hh"
#include "test_solid_mechanics_model_fixture.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

namespace {

template <typename type_>
class TestSMMFixtureDynamics : public TestSMMFixture<type_> {

public:
  static constexpr ElementType type = type_::value;

  void SetUp() override {
    if (this->type == _pentahedron_6 || this->type == _pentahedron_15)
      return;

    std::cout << "testing type " << this->type << std::endl;

    TestSMMFixture<type_>::SetUp();
  }

  std::string makeMeshName() override {
    std::stringstream element_type;
    element_type << type;
    SCOPED_TRACE(element_type.str().c_str());
    return std::string("bar") + element_type.str() + ".msh";
  }
};

constexpr Real A = 1e-1;
constexpr Real E = 1.;
constexpr Real poisson = .3;
constexpr Real lambda = E * poisson / (1. + poisson) / (1. - 2. * poisson);
constexpr Real mu = E / 2 / (1. + poisson);
constexpr Real rho = 1.;
const Real cp = std::sqrt((lambda + 2 * mu) / rho);
const Real cs = std::sqrt(mu / rho);
const Real c = std::sqrt(E / rho);

const Vector<Real> k = {.5, 0., 0.};
const Vector<Real> psi1 = {0., 0., 1.};
const Vector<Real> psi2 = {0., 1., 0.};
const Real knorm = k.norm();

template <UInt spatial_dimension> auto solution_disp = 1;
template <UInt spatial_dimension> auto solution_vel = 1;

template <>
auto solution_disp<1> =
    [](Vector<Real> & disp, const Vector<Real> & coord, Real current_time) {
      const auto & x = coord(_x);
      const Real omega = c * k[0];
      disp(_x) = A * cos(k[0] * x - omega * current_time);
    };

template <>
auto solution_vel<1> =
    [](Vector<Real> & vel, const Vector<Real> & coord, Real current_time) {
      const auto & x = coord(_x);
      const Real omega = c * k[0];
      vel(_x) = omega * A * sin(k[0] * x - omega * current_time);
    };

template <>
auto solution_disp<2> =
    [](Vector<Real> & disp, const Vector<Real> & coord, Real current_time) {

      const auto & X = coord;
      Vector<Real> kshear = {k[1], k[0]};
      Vector<Real> kpush = {k[0], k[1]};

      const Real omega_p = knorm * cp;
      const Real omega_s = knorm * cs;

      Real phase_p = X.dot(kpush) - omega_p * current_time;
      Real phase_s = X.dot(kpush) - omega_s * current_time;

      disp = A * (kpush * cos(phase_p) + kshear * cos(phase_s));
    };

template <>
auto solution_vel<2> = [](Vector<Real> & vel, const Vector<Real> & coord,
                          Real current_time) {

  const auto & X = coord;
  Vector<Real> kshear = {k[1], k[0]};
  Vector<Real> kpush = {k[0], k[1]};

  const Real omega_p = knorm * cp;
  const Real omega_s = knorm * cs;

  Real phase_p = X.dot(kpush) - omega_p * current_time;
  Real phase_s = X.dot(kpush) - omega_s * current_time;

  vel = A * (kpush * omega_p * cos(phase_p) + kshear * omega_s * cos(phase_s));

};

template <>
auto solution_disp<3> =
    [](Vector<Real> & disp, const Vector<Real> & coord, Real current_time) {

      const auto & X = coord;
      Vector<Real> kpush = k;
      Vector<Real> kshear1(3);
      Vector<Real> kshear2(3);
      kshear1.crossProduct(k, psi1);
      kshear2.crossProduct(k, psi2);

      const Real omega_p = knorm * cp;
      const Real omega_s = knorm * cs;

      Real phase_p = X.dot(kpush) - omega_p * current_time;
      Real phase_s = X.dot(kpush) - omega_s * current_time;

      disp = A * (kpush * cos(phase_p) + kshear1 * cos(phase_s) +
                  kshear2 * cos(phase_s));

    };

template <>
auto solution_vel<3> = [](Vector<Real> & vel, const Vector<Real> & coord,
                          Real current_time) {

  const auto & X = coord;
  Vector<Real> kpush = k;
  Vector<Real> kshear1(3);
  Vector<Real> kshear2(3);
  kshear1.crossProduct(k, psi1);
  kshear2.crossProduct(k, psi2);

  const Real omega_p = knorm * cp;
  const Real omega_s = knorm * cs;

  Real phase_p = X.dot(kpush) - omega_p * current_time;
  Real phase_s = X.dot(kpush) - omega_s * current_time;

  vel = A * (kpush * omega_p * cos(phase_p) + kshear1 * omega_s * cos(phase_s) +
             kshear2 * omega_s * cos(phase_s));
};

template <ElementType _type>
class SolutionFunctor : public BC::Dirichlet::DirichletFunctor {
public:
  SolutionFunctor(Real current_time, SolidMechanicsModel & model)
      : current_time(current_time), model(model) {}

public:
  // static constexpr UInt dim = DimensionHelper<_type>::dim;
  static constexpr UInt dim = ElementClass<_type>::getSpatialDimension();

  inline void operator()(UInt node, Vector<bool> & flags, Vector<Real> & primal,
                         const Vector<Real> & coord) const {

    flags(0) = true;
    auto & vel = model.getVelocity();
    auto itVel = vel.begin(model.getSpatialDimension());
    Vector<Real> v = itVel[node];
    solution_disp<dim>(primal, coord, current_time);
    solution_vel<dim>(v, coord, current_time);
  }

private:
  Real current_time;
  SolidMechanicsModel & model;
};

template <ElementType _type, typename AM>
void test_body(SolidMechanicsModel & model, AM analysis_method) {

  // constexpr UInt dim = DimensionHelper<_type>::dim;
  static constexpr UInt dim = ElementClass<_type>::getSpatialDimension();

  getStaticParser().parse("test_solid_mechanics_model_"
                          "dynamics_material.dat");

  model.initFull(SolidMechanicsModelOptions(analysis_method));

  bool dump_paraview = false;

  if (dump_paraview) {
    std::stringstream base_name;
    base_name << "bar" << analysis_method << _type;
    model.setBaseName(base_name.str());
    model.addDumpFieldVector("displacement");
    model.addDumpField("mass");
    model.addDumpField("velocity");
    model.addDumpField("acceleration");
    model.addDumpFieldVector("external_force");
    model.addDumpFieldVector("internal_force");
    model.addDumpField("stress");
    model.addDumpField("strain");
  }

  Real time_step = model.getStableTimeStep() / 10.;
  model.setTimeStep(time_step);
  std::cout << "timestep: " << time_step << std::endl;

  UInt nb_nodes = model.getMesh().getNbNodes();
  UInt spatial_dimension = model.getSpatialDimension();

  auto & nodes = model.getMesh().getNodes();
  auto & disp = model.getDisplacement();
  auto & vel = model.getVelocity();

  Array<Real> disp_solution(nb_nodes, spatial_dimension);

  Real current_time = 0;

  auto itNodes = nodes.begin(spatial_dimension);
  auto itDisp = disp.begin(spatial_dimension);
  auto itVel = vel.begin(spatial_dimension);
  for (UInt n = 0; n < nb_nodes; ++n, ++itNodes, ++itDisp, ++itVel) {
    solution_disp<dim>(*itDisp, *itNodes, current_time);
    solution_vel<dim>(*itVel, *itNodes, current_time);
  }

  if (dump_paraview)
    model.dump();

  /// boundary conditions
  model.applyBC(SolutionFunctor<_type>(current_time, model), "BC");

  Real max_error = 0.;
  Real wave_velocity = 1.; // sqrt(E/rho) = sqrt(1/1) = 1
  Real simulation_time = 5 / wave_velocity;

  UInt max_steps = simulation_time / time_step; // 100
  UInt ndump = 50;
  UInt dump_freq = max_steps / ndump;

  std::cout << "max_steps: " << max_steps << std::endl;

  for (UInt s = 0; s < max_steps; ++s, current_time += time_step) {

    if (s % dump_freq == 0 && dump_paraview)
      model.dump();

    /// boundary conditions
    model.applyBC(SolutionFunctor<_type>(current_time, model), "BC");

    // compute the disp solution
    auto itDispSolution = disp_solution.begin(spatial_dimension);
    itNodes = nodes.begin(spatial_dimension);
    for (UInt n = 0; n < nb_nodes; ++n, ++itNodes, ++itDispSolution) {
      solution_disp<dim>(*itDispSolution, *itNodes, current_time);
    }
    // compute the error solution
    itDispSolution = disp_solution.begin(spatial_dimension);
    itDisp = disp.begin(spatial_dimension);
    Real disp_error = 0.;
    for (UInt n = 0; n < nb_nodes; ++n, ++itDispSolution, ++itDisp) {
      auto diff = *itDispSolution - *itDisp;

      for (UInt i = 0; i < spatial_dimension; ++i) {
        disp_error += diff(i) * diff(i);
      }
    }
    disp_error = sqrt(disp_error) / nb_nodes;
    max_error = std::max(disp_error, max_error);
    ASSERT_NEAR(disp_error, 0., 1e-3);
    model.solveStep();
  }
  std::cout << "max error: " << max_error << std::endl;
}

TYPED_TEST_CASE(TestSMMFixtureDynamics, types);

#ifdef AKANTU_IMPLICIT
TYPED_TEST(TestSMMFixture, DynamicsImplicit) {
  if (this->type != _pentahedron_6 && this->type != _pentahedron_15)
    test_body<this->type>(*(this->model), _implicit_dynamic);
}
#endif

TYPED_TEST(TestSMMFixtureDynamics, DynamicsExplicit) {
  if (this->type != _pentahedron_6 && this->type != _pentahedron_15)
    test_body<this->type>(*(this->model), _explicit_lumped_mass);
}
}
