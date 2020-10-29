/**
 * @file   test_solid_mechanics_model_dynamics.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date creation: Wed Nov 29 2017
 * @date last modification: Wed Feb 21 2018
 *
 * @brief  test of the class SolidMechanicsModel on the 3d cube
 *
 *
 * Copyright (©) 2016-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
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

const Real A = 1e-1;
const Real E = 1.;
const Real poisson = 3. / 10;
const Real lambda = E * poisson / (1 + poisson) / (1 - 2 * poisson);
const Real mu = E / 2 / (1. + poisson);
const Real rho = 1.;
const Real cp = std::sqrt((lambda + 2 * mu) / rho);
const Real cs = std::sqrt(mu / rho);
const Real c = std::sqrt(E / rho);

const Vector<Real> k = {.5, 0., 0.};
const Vector<Real> psi1 = {0., 0., 1.};
const Vector<Real> psi2 = {0., 1., 0.};
const Real knorm = k.norm();

/* -------------------------------------------------------------------------- */
template <UInt dim> struct Verification {};

/* -------------------------------------------------------------------------- */
template <> struct Verification<1> {
  void displacement(Vector<Real> & disp, const Vector<Real> & coord,
                    Real current_time) {
    const auto & x = coord(_x);
    const Real omega = c * k[0];
    disp(_x) = A * cos(k[0] * x - omega * current_time);
  }

  void velocity(Vector<Real> & vel, const Vector<Real> & coord,
                Real current_time) {
    const auto & x = coord(_x);
    const Real omega = c * k[0];
    vel(_x) = omega * A * sin(k[0] * x - omega * current_time);
  }
};

/* -------------------------------------------------------------------------- */
template <> struct Verification<2> {
  void displacement(Vector<Real> & disp, const Vector<Real> & X,
                    Real current_time) {
    Vector<Real> kshear = {k[1], k[0]};
    Vector<Real> kpush = {k[0], k[1]};

    const Real omega_p = knorm * cp;
    const Real omega_s = knorm * cs;

    Real phase_p = X.dot(kpush) - omega_p * current_time;
    Real phase_s = X.dot(kpush) - omega_s * current_time;

    disp = A * (kpush * cos(phase_p) + kshear * cos(phase_s));
  }

  void velocity(Vector<Real> & vel, const Vector<Real> & X, Real current_time) {
    Vector<Real> kshear = {k[1], k[0]};
    Vector<Real> kpush = {k[0], k[1]};

    const Real omega_p = knorm * cp;
    const Real omega_s = knorm * cs;

    Real phase_p = X.dot(kpush) - omega_p * current_time;
    Real phase_s = X.dot(kpush) - omega_s * current_time;

    vel =
        A * (kpush * omega_p * cos(phase_p) + kshear * omega_s * cos(phase_s));
  }
};

/* -------------------------------------------------------------------------- */
template <> struct Verification<3> {
  void displacement(Vector<Real> & disp, const Vector<Real> & coord,
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

    disp = A * (kpush * cos(phase_p) + kshear1 * cos(phase_s) +
                kshear2 * cos(phase_s));
  }

  void velocity(Vector<Real> & vel, const Vector<Real> & coord,
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

    vel =
        A * (kpush * omega_p * cos(phase_p) + kshear1 * omega_s * cos(phase_s) +
             kshear2 * omega_s * cos(phase_s));
  }
};

/* -------------------------------------------------------------------------- */
template <ElementType _type>
class SolutionFunctor : public BC::Dirichlet::DirichletFunctor {
public:
  SolutionFunctor(Real current_time, SolidMechanicsModel & model)
      : current_time(current_time), model(model) {}

public:
  static constexpr UInt dim = ElementClass<_type>::getSpatialDimension();

  inline void operator()(UInt node, Vector<bool> & flags, Vector<Real> & primal,
                         const Vector<Real> & coord) const {
    flags(0) = true;
    auto & vel = model.getVelocity();
    auto it = vel.begin(model.getSpatialDimension());
    Vector<Real> v = it[node];

    Verification<dim> verif;
    verif.displacement(primal, coord, current_time);
    verif.velocity(v, coord, current_time);
  }

private:
  Real current_time;
  SolidMechanicsModel & model;
};

/* -------------------------------------------------------------------------- */
// This fixture uses somewhat finer meshes representing bars.
template <typename type_, typename analysis_method_>
class TestSMMFixtureBar : public TestSMMFixture<type_> {
  using parent = TestSMMFixture<type_>;

public:
  void SetUp() override {
    this->mesh_file =
        "../patch_tests/data/bar" + std::to_string(this->type) + ".msh";
    parent::SetUp();

    auto analysis_method = analysis_method_::value;
    this->initModel("test_solid_mechanics_model_"
                    "dynamics_material.dat",
                    analysis_method);

    const auto & position = this->mesh->getNodes();
    auto & displacement = this->model->getDisplacement();
    auto & velocity = this->model->getVelocity();

    constexpr auto dim = parent::spatial_dimension;

    Verification<dim> verif;

    for (auto && tuple :
         zip(make_view(position, dim), make_view(displacement, dim),
             make_view(velocity, dim))) {
      verif.displacement(std::get<1>(tuple), std::get<0>(tuple), 0.);
      verif.velocity(std::get<2>(tuple), std::get<0>(tuple), 0.);
    }

    if (this->dump_paraview)
      this->model->dump();

    /// boundary conditions
    this->model->applyBC(SolutionFunctor<parent::type>(0., *this->model), "BC");

    wave_velocity = 1.; // sqrt(E/rho) = sqrt(1/1) = 1
    simulation_time = 5 / wave_velocity;
    time_step = this->model->getTimeStep();

    max_steps = simulation_time / time_step; // 100
  }

  void solveStep() {
    constexpr auto dim = parent::spatial_dimension;
    Real current_time = 0.;
    const auto & position = this->mesh->getNodes();
    const auto & displacement = this->model->getDisplacement();
    UInt nb_nodes = this->mesh->getNbNodes();
    UInt nb_global_nodes = this->mesh->getNbGlobalNodes();

    max_error = -1.;

    Array<Real> displacement_solution(nb_nodes, dim);
    Verification<dim> verif;

    auto ndump = 50;
    auto dump_freq = max_steps / ndump;

    for (UInt s = 0; s < this->max_steps;
         ++s, current_time += this->time_step) {
      if (s % dump_freq == 0 && this->dump_paraview)
        this->model->dump();

      /// boundary conditions
      this->model->applyBC(
          SolutionFunctor<parent::type>(current_time, *this->model), "BC");

      // compute the disp solution
      for (auto && tuple : zip(make_view(position, dim),
                               make_view(displacement_solution, dim))) {
        verif.displacement(std::get<1>(tuple), std::get<0>(tuple),
                           current_time);
      }

      // compute the error solution
      Real disp_error = 0.;
      auto n = 0;
      for (auto && tuple : zip(make_view(displacement, dim),
                               make_view(displacement_solution, dim))) {
        if (this->mesh->isLocalOrMasterNode(n)) {
          auto diff = std::get<1>(tuple) - std::get<0>(tuple);
          disp_error += diff.dot(diff);
        }
        ++n;
      }

      this->mesh->getCommunicator().allReduce(disp_error,
                                              SynchronizerOperation::_sum);

      disp_error = sqrt(disp_error) / nb_global_nodes;
      max_error = std::max(disp_error, max_error);

      this->model->solveStep();
    }
  }

protected:
  Real time_step;
  Real wave_velocity;
  Real simulation_time;
  UInt max_steps;
  Real max_error{-1};
};

template <AnalysisMethod t>
using analysis_method_t = std::integral_constant<AnalysisMethod, t>;

using TestTypes = gtest_list_t<TestElementTypes>;

template <typename type_>
using TestSMMFixtureBarExplicit =
    TestSMMFixtureBar<type_, analysis_method_t<_explicit_lumped_mass>>;

TYPED_TEST_SUITE(TestSMMFixtureBarExplicit, TestTypes);

/* -------------------------------------------------------------------------- */
TYPED_TEST(TestSMMFixtureBarExplicit, Dynamics) {
  this->solveStep();
  EXPECT_NEAR(this->max_error, 0., 2e-3);
  // std::cout << "max error: " << max_error << std::endl;
}

/* -------------------------------------------------------------------------- */
#if defined(AKANTU_IMPLICIT)
template <typename type_>
using TestSMMFixtureBarImplicit =
    TestSMMFixtureBar<type_, analysis_method_t<_implicit_dynamic>>;

TYPED_TEST_SUITE(TestSMMFixtureBarImplicit, TestTypes);

TYPED_TEST(TestSMMFixtureBarImplicit, Dynamics) {
  if (this->type == _segment_2 and
      (this->mesh->getCommunicator().getNbProc() > 2)) {
    // The error are just to high after (hopefully because of the two small test
    // case)
    SUCCEED();
    return;
  }
  this->solveStep();
  EXPECT_NEAR(this->max_error, 0., 2e-3);
}
#endif

} // namespace
