/**
 * @file   test_lumped_mass.cc
 *
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Dec 05 2017
 * @date last modification: Tue Jan 30 2018
 *
 * @brief  test the lumping of the mass matrix
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
#include "solid_mechanics_model.hh"
#include "test_gtest_utils.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
#include <tuple>
/* -------------------------------------------------------------------------- */

using namespace akantu;

template <typename tuple_>
class TestLumpedMassesFixture : public ::testing::Test {
public:
  static constexpr ElementType type = tuple_::value;
  static constexpr size_t dim = ElementClass<type>::getSpatialDimension();

  void SetUp() override {
    debug::setDebugLevel(dblError);
    getStaticParser().parse("material_lumped.dat");

    std::stringstream element_type;
    element_type << type;

    mesh = std::make_unique<Mesh>(dim);
    mesh->read(element_type.str() + ".msh");

    SCOPED_TRACE(element_type.str().c_str());

    model = std::make_unique<SolidMechanicsModel>(*mesh);

    model->initFull(_analysis_method = _explicit_lumped_mass);
  }

  void TearDown() override {
    model.reset(nullptr);
    mesh.reset(nullptr);
  }

protected:
  std::unique_ptr<Mesh> mesh;
  std::unique_ptr<SolidMechanicsModel> model;
};

template <typename T> constexpr ElementType TestLumpedMassesFixture<T>::type;
template <typename T> constexpr size_t TestLumpedMassesFixture<T>::dim;

using mass_types = gtest_list_t<TestElementTypes>;

TYPED_TEST_SUITE(TestLumpedMassesFixture, mass_types);

TYPED_TEST(TestLumpedMassesFixture, TestLumpedMass) {
  this->model->assembleMassLumped();

  auto rho = this->model->getMaterial(0).getRho();
  auto & fem = this->model->getFEEngine();
  auto nb_element = this->mesh->getNbElement(this->type);
  auto nb_quadrature_points =
      fem.getNbIntegrationPoints(this->type) * nb_element;

  Array<Real> rho_on_quad(nb_quadrature_points, 1, rho, "rho_on_quad");
  auto mass = fem.integrate(rho_on_quad, this->type);
  const auto & masses = this->model->getMass();

  Vector<Real> sum(this->dim, 0.);
  for (auto & mass : make_view(masses, this->dim)) {
    sum += mass;
  }

  for (UInt s = 0; s < sum.size(); ++s)
    EXPECT_NEAR(0., (mass - sum[s]) / mass, 2e-15);
}
