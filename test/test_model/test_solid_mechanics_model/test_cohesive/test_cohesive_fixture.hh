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

template <typename param_> class TestSMMCFixture : public ::testing::Test {
public:
  static constexpr ElementType cohesive_type =
      std::tuple_element_t<0, param_>::value;
  static constexpr ElementType type_1 = std::tuple_element_t<1, param_>::value;
  static constexpr ElementType type_2 = std::tuple_element_t<2, param_>::value;

  static constexpr size_t dim =
      ElementClass<cohesive_type>::getSpatialDimension();

  void SetUp() override {
    normal = Vector<Real>(this->dim, 0.);
    if (dim == 1)
      normal(_x) = 1.;
    else
      normal(_y) = 1.;

    mesh = std::make_unique<Mesh>(this->dim);
    if (Communicator::getStaticCommunicator().whoAmI() == 0) {
      EXPECT_NO_THROW({ mesh->read(this->mesh_name); });
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

    //    auto stable_time_step = this->model->getStableTimeStep();
    this->model->setTimeStep(4e-6);

    //  std::cout << stable_time_step *0.0 << std::endl;
    if (dim == 1) {
      surface = 1;
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
  }

  void setInitialCondition(const Vector<Real> & direction, Real strain) {
    auto lower = this->mesh->getLowerBounds().dot(normal);
    auto upper = this->mesh->getUpperBounds().dot(normal);

    Real L = upper - lower;

    for (auto && data :
         zip(make_view(this->mesh->getNodes(), this->dim),
             make_view(this->model->getDisplacement(), this->dim))) {
      const auto & pos = std::get<0>(data);
      auto & disp = std::get<1>(data);

      auto x = pos.dot(normal) - (upper + lower) / 2.;
      disp = direction * (x * 2. * strain / L);
    }
  }

#define debug 1

  void steps(const Vector<Real> & displacement_max) {
#if debug
    this->model->addDumpFieldVector("displacement");
    this->model->addDumpFieldVector("velocity");
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
    auto inc_load = BC::Dirichlet::Increment(displacement_max / Real(nb_steps));
    auto inc_fix =
        BC::Dirichlet::Increment(-1. / Real(nb_steps) * displacement_max);

    for (auto _[[gnu::unused]] : arange(nb_steps)) {
      this->model->applyBC(inc_load, "loading");
      this->model->applyBC(inc_fix, "fixed");
      if (this->is_extrinsic)
        this->model->checkCohesiveStress();

      this->model->solveStep();
#if debug
      this->model->dump();
      this->model->dump("cohesive elements");
#endif
    }
  }
  void checkInsertion() {
    auto nb_cohesive_element = this->mesh->getNbElement(cohesive_type);
    mesh->getCommunicator().allReduce(nb_cohesive_element,
                                      SynchronizerOperation::_sum);

    auto facet_type = this->mesh->getFacetType(this->cohesive_type);
    const auto & group =
        this->mesh->getElementGroup("insertion").getElements(facet_type);
    auto group_size = group.size();

    mesh->getCommunicator().allReduce(group_size, SynchronizerOperation::_sum);

    EXPECT_EQ(nb_cohesive_element, group_size);
  }

  void checkDissipated(Real expected_density) {
    Real edis = this->model->getEnergy("dissipated");

    EXPECT_NEAR(this->surface * expected_density, edis, 4e-1);
  }

  void testModeI() {
    Vector<Real> direction(this->dim, 0.);
    if (dim == 1)
      direction(_x) = 1.;
    else
      direction(_y) = 1.;

    //  EXPECT_NO_THROW(this->createModel());
    this->createModel();

    if (this->dim > 1)
      this->model->applyBC(BC::Dirichlet::FlagOnly(_x), "sides");
    if (this->dim > 2)
      this->model->applyBC(BC::Dirichlet::FlagOnly(_z), "sides");

    auto & mat_co = this->model->getMaterial("insertion");
    Real sigma_c = mat_co.get("sigma_c");
    Real G_c = mat_co.get("G_c");

    auto & mat_el = this->model->getMaterial("body");
    Real E = mat_el.get("E");
    Real nu = mat_el.get("nu");

    auto delta_c = 2. * G_c / sigma_c;
    Real strain = sigma_c / E;

    if (dim == 1)
      strain *= .9999;
    else
      strain *= .935; // there must be an error in my computations

    if (this->dim == 2)
      strain *= (1. - nu) * (1. + nu);

    auto max_travel = .5 * delta_c;
    this->setInitialCondition(direction, strain);
    this->steps(direction * max_travel);
  }

  void testModeII() {
    Vector<Real> direction(this->dim, 0.);
    direction(_x) = 1.;

    EXPECT_NO_THROW(this->createModel());

    if (this->dim > 1)
      this->model->applyBC(BC::Dirichlet::FlagOnly(_y), "sides");
    if (this->dim > 2)
      this->model->applyBC(BC::Dirichlet::FlagOnly(_z), "sides");

    auto & mat_co = this->model->getMaterial("insertion");
    Real sigma_c = mat_co.get("sigma_c");
    Real G_c = mat_co.get("G_c");
    Real beta = mat_co.get("beta");

    auto & mat_el = this->model->getMaterial("body");
    Real E = mat_el.get("E");
    Real nu = mat_el.get("nu");

    auto L = this->mesh->getUpperBounds().dot(direction) -
             this->mesh->getLowerBounds().dot(direction);

    auto delta_c = 2. * G_c / sigma_c;
    Real strain = .99999 * L * beta * beta * sigma_c / E;
    if (this->dim > 1) {
      strain *= (1. + nu);
    }
    auto max_travel = 1.2 * delta_c;
    this->setInitialCondition(direction, strain);
    this->steps(direction * max_travel);
  }

protected:
  std::unique_ptr<Mesh> mesh;
  std::unique_ptr<SolidMechanicsModelCohesive> model;

  std::string mesh_name{aka::to_string(cohesive_type) + aka::to_string(type_1) +
                        (type_1 == type_2 ? "" : aka::to_string(type_2)) +
                        ".msh"};

  bool is_extrinsic;
  AnalysisMethod analysis_method;

  Vector<Real> normal;
  Real surface{0};
  UInt nb_steps{300};
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

using types = gtest_list_t<std::tuple<
    std::tuple<element_type_t<_cohesive_1d_2>, element_type_t<_segment_2>,
               element_type_t<_segment_2>>,
    std::tuple<element_type_t<_cohesive_2d_4>, element_type_t<_triangle_3>,
               element_type_t<_triangle_3>>,
    std::tuple<element_type_t<_cohesive_2d_4>, element_type_t<_quadrangle_4>,
               element_type_t<_quadrangle_4>>,
    std::tuple<element_type_t<_cohesive_2d_4>, element_type_t<_triangle_3>,
               element_type_t<_quadrangle_4>>,
    std::tuple<element_type_t<_cohesive_2d_6>, element_type_t<_triangle_6>,
               element_type_t<_triangle_6>>,
    std::tuple<element_type_t<_cohesive_2d_6>, element_type_t<_quadrangle_8>,
               element_type_t<_quadrangle_8>>,
    std::tuple<element_type_t<_cohesive_2d_6>, element_type_t<_triangle_6>,
               element_type_t<_quadrangle_8>>,
    std::tuple<element_type_t<_cohesive_3d_6>, element_type_t<_tetrahedron_4>,
               element_type_t<_tetrahedron_4>>,
    std::tuple<element_type_t<_cohesive_3d_12>, element_type_t<_tetrahedron_10>,
               element_type_t<_tetrahedron_10>>,
    std::tuple<element_type_t<_cohesive_3d_8>, element_type_t<_hexahedron_8>,
               element_type_t<_hexahedron_8>>,
    std::tuple<element_type_t<_cohesive_3d_16>, element_type_t<_hexahedron_20>,
               element_type_t<_hexahedron_20>>>>;

TYPED_TEST_CASE(TestSMMCFixture, types);

#endif /* __AKANTU_TEST_COHESIVE_FIXTURE_HH__ */
