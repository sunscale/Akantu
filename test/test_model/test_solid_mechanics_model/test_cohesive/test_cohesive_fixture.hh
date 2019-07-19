/**
 * @file   test_cohesive_fixture.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Jan 10 2018
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Coehsive element test fixture
 *
 * @section LICENSE
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
#include "communicator.hh"
#include "solid_mechanics_model_cohesive.hh"
#include "test_gtest_utils.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
#include <vector>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_TEST_COHESIVE_FIXTURE_HH__
#define __AKANTU_TEST_COHESIVE_FIXTURE_HH__

using namespace akantu;

template <::akantu::AnalysisMethod t>
using analysis_method_t = std::integral_constant<::akantu::AnalysisMethod, t>;

class StrainIncrement : public BC::Functor {
public:
  StrainIncrement(const Matrix<Real> & strain, BC::Axis dir)
      : strain_inc(strain), dir(dir) {}

  void operator()(UInt /*node*/, Vector<bool> & flags, Vector<Real> & primal,
                  const Vector<Real> & coord) const {
    if (std::abs(coord(dir)) < 1e-8) {
      return;
    }

    flags.set(true);
    primal += strain_inc * coord;
  }

  static const BC::Functor::Type type = BC::Functor::_dirichlet;

private:
  Matrix<Real> strain_inc;
  BC::Axis dir;
};

template <typename param_> class TestSMMCFixture : public ::testing::Test {
public:
  static constexpr ElementType cohesive_type =
      std::tuple_element_t<0, param_>::value;
  static constexpr ElementType type_1 = std::tuple_element_t<1, param_>::value;
  static constexpr ElementType type_2 = std::tuple_element_t<2, param_>::value;

  static constexpr size_t dim =
      ElementClass<cohesive_type>::getSpatialDimension();

  void SetUp() override {
    mesh = std::make_unique<Mesh>(this->dim);
    if (Communicator::getStaticCommunicator().whoAmI() == 0) {
      ASSERT_NO_THROW({ mesh->read(this->mesh_name); });
    }
    mesh->distribute();
  }

  void TearDown() override {
    model.reset(nullptr);
    mesh.reset(nullptr);
  }

  void createModel() {
    model = std::make_unique<SolidMechanicsModelCohesive>(*mesh);
    model->initFull(_analysis_method = this->analysis_method,
                    _is_extrinsic = this->is_extrinsic);

    auto time_step = this->model->getStableTimeStep() * 0.01;
    this->model->setTimeStep(time_step);

    if (dim == 1) {
      surface = 1;
      group_size = 1;
      return;
    }

    auto facet_type = mesh->getFacetType(this->cohesive_type);

    auto & fe_engine = model->getFEEngineBoundary();
    auto & group = mesh->getElementGroup("insertion").getElements(facet_type);
    Array<Real> ones(fe_engine.getNbIntegrationPoints(facet_type) *
                     group.size());
    ones.set(1.);

    surface = fe_engine.integrate(ones, facet_type, _not_ghost, group);
    mesh->getCommunicator().allReduce(surface, SynchronizerOperation::_sum);

    group_size = group.size();

    mesh->getCommunicator().allReduce(group_size, SynchronizerOperation::_sum);

#define debug_ 0

#if debug_
    this->model->addDumpFieldVector("displacement");
    this->model->addDumpFieldVector("velocity");
    this->model->addDumpFieldVector("internal_force");
    this->model->addDumpFieldVector("external_force");
    this->model->addDumpField("blocked_dofs");
    this->model->addDumpField("stress");
    this->model->addDumpField("strain");
    this->model->assembleInternalForces();
    this->model->setBaseNameToDumper("cohesive elements", "cohesive_elements");
    this->model->addDumpFieldVectorToDumper("cohesive elements",
                                            "displacement");
    this->model->addDumpFieldToDumper("cohesive elements", "damage");
    this->model->addDumpFieldToDumper("cohesive elements", "tractions");
    this->model->addDumpFieldToDumper("cohesive elements", "opening");
    this->model->dump();
    this->model->dump("cohesive elements");
#endif
  }

  void setInitialCondition(const Matrix<Real> & strain) {
    for (auto && data :
         zip(make_view(this->mesh->getNodes(), this->dim),
             make_view(this->model->getDisplacement(), this->dim))) {
      const auto & pos = std::get<0>(data);
      auto & disp = std::get<1>(data);
      disp = strain * pos;
    }
  }

  bool checkDamaged() {
    UInt nb_damaged = 0;
    auto & damage =
        model->getMaterial("insertion").getArray<Real>("damage", cohesive_type);
    for (auto d : damage) {
      if (d >= .99)
        ++nb_damaged;
    }

    return (nb_damaged == group_size);
  }

  void steps(const Matrix<Real> & strain) {
    StrainIncrement functor((1. / 300) * strain, this->dim == 1 ? _x : _y);

    for (auto _ [[gnu::unused]] : arange(nb_steps)) {
      this->model->applyBC(functor, "loading");
      this->model->applyBC(functor, "fixed");
      if (this->is_extrinsic)
        this->model->checkCohesiveStress();

      this->model->solveStep();
#if debug_
      this->model->dump();
      this->model->dump("cohesive elements");
#endif
    }
  }

  void checkInsertion() {
    auto nb_cohesive_element = this->mesh->getNbElement(cohesive_type);
    mesh->getCommunicator().allReduce(nb_cohesive_element,
                                      SynchronizerOperation::_sum);

    EXPECT_EQ(nb_cohesive_element, group_size);
  }

  void checkDissipated(Real expected_density) {
    Real edis = this->model->getEnergy("dissipated");

    EXPECT_NEAR(this->surface * expected_density, edis, 5e-1);
  }

  void testModeI() {
    this->createModel();

    auto & mat_el = this->model->getMaterial("body");

    auto speed = mat_el.getPushWaveSpeed(Element());
    auto direction = _y;
    if (dim == 1)
      direction = _x;
    auto length =
        mesh->getUpperBounds()(direction) - mesh->getLowerBounds()(direction);
    nb_steps = length / 2. / speed / model->getTimeStep();

    SCOPED_TRACE(std::to_string(this->dim) + "D - " + std::to_string(type_1) +
                 ":" + std::to_string(type_2));

    auto & mat_co = this->model->getMaterial("insertion");
    Real sigma_c = mat_co.get("sigma_c");

    Real E = mat_el.get("E");
    Real nu = mat_el.get("nu");

    Matrix<Real> strain;
    if (dim == 1) {
      strain = {{1.}};
    } else if (dim == 2) {
      strain = {{-nu, 0.}, {0., 1. - nu}};
      strain *= (1. + nu);
    } else if (dim == 3) {
      strain = {{-nu, 0., 0.}, {0., 1., 0.}, {0., 0., -nu}};
    }

    strain *= sigma_c / E;

    this->setInitialCondition((1 - 1e-5) * strain);
    this->steps(1e-2 * strain);
  }

  void testModeII() {
    this->createModel();
    auto & mat_el = this->model->getMaterial("body");
    Real speed;
    try {
      speed =
          mat_el.getShearWaveSpeed(Element()); // the slowest speed if exists
    } catch (...) {
      speed = mat_el.getPushWaveSpeed(Element());
    }

    auto direction = _y;
    if (dim == 1)
      direction = _x;
    auto length =
        mesh->getUpperBounds()(direction) - mesh->getLowerBounds()(direction);
    nb_steps = 2 * length / 2. / speed / model->getTimeStep();

    SCOPED_TRACE(std::to_string(this->dim) + "D - " + std::to_string(type_1) +
                 ":" + std::to_string(type_2));

    if (this->dim > 1)
      this->model->applyBC(BC::Dirichlet::FlagOnly(_y), "sides");
    if (this->dim > 2)
      this->model->applyBC(BC::Dirichlet::FlagOnly(_z), "sides");

    auto & mat_co = this->model->getMaterial("insertion");
    Real sigma_c = mat_co.get("sigma_c");
    Real beta = mat_co.get("beta");
    // Real G_c = mat_co.get("G_c");

    Real E = mat_el.get("E");
    Real nu = mat_el.get("nu");

    Matrix<Real> strain;
    if (dim == 1) {
      strain = {{1.}};
    } else if (dim == 2) {
      strain = {{0., 1.}, {0., 0.}};
      strain *= (1. + nu);
    } else if (dim == 3) {
      strain = {{0., 1., 0.}, {0., 0., 0.}, {0., 0., 0.}};
      strain *= (1. + nu);
    }
    strain *= 2 * beta * beta * sigma_c / E;

    // nb_steps *= 5;

    this->setInitialCondition((1. - 1e-5) * strain);
    this->steps(0.005 * strain);
  }

protected:
  std::unique_ptr<Mesh> mesh;
  std::unique_ptr<SolidMechanicsModelCohesive> model;

  std::string mesh_name{std::to_string(cohesive_type) + std::to_string(type_1) +
                        (type_1 == type_2 ? "" : std::to_string(type_2)) +
                        ".msh"};

  bool is_extrinsic;
  AnalysisMethod analysis_method;

  Real surface{0};
  UInt nb_steps{1000};
  UInt group_size{10000};
};

/* -------------------------------------------------------------------------- */
template <typename param_>
constexpr ElementType TestSMMCFixture<param_>::cohesive_type;
template <typename param_>
constexpr ElementType TestSMMCFixture<param_>::type_1;
template <typename param_>
constexpr ElementType TestSMMCFixture<param_>::type_2;

template <typename param_> constexpr size_t TestSMMCFixture<param_>::dim;
/* -------------------------------------------------------------------------- */

using IsExtrinsicTypes = std::tuple<std::true_type, std::false_type>;
using AnalysisMethodTypes =
    std::tuple<analysis_method_t<_explicit_lumped_mass>>;

using coh_types = gtest_list_t<std::tuple<
    std::tuple<_element_type_cohesive_1d_2, _element_type_segment_2,
               _element_type_segment_2>,
    std::tuple<_element_type_cohesive_2d_4, _element_type_triangle_3,
               _element_type_triangle_3>,
    std::tuple<_element_type_cohesive_2d_4, _element_type_quadrangle_4,
               _element_type_quadrangle_4>,
    std::tuple<_element_type_cohesive_2d_4, _element_type_triangle_3,
               _element_type_quadrangle_4>,
    std::tuple<_element_type_cohesive_2d_6, _element_type_triangle_6,
               _element_type_triangle_6>,
    std::tuple<_element_type_cohesive_2d_6, _element_type_quadrangle_8,
               _element_type_quadrangle_8>,
    std::tuple<_element_type_cohesive_2d_6, _element_type_triangle_6,
               _element_type_quadrangle_8>,
    std::tuple<_element_type_cohesive_3d_6, _element_type_tetrahedron_4,
               _element_type_tetrahedron_4>,
    std::tuple<_element_type_cohesive_3d_12, _element_type_tetrahedron_10,
               _element_type_tetrahedron_10> /*,
    std::tuple<_element_type_cohesive_3d_8, _element_type_hexahedron_8,
               _element_type_hexahedron_8>,
    std::tuple<_element_type_cohesive_3d_16, _element_type_hexahedron_20,
               _element_type_hexahedron_20>*/>>;

TYPED_TEST_SUITE(TestSMMCFixture, coh_types);

#endif /* __AKANTU_TEST_COHESIVE_FIXTURE_HH__ */
