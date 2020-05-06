/**
 * @file   test_damage_materials.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date creation: Fri Nov 17 2017
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Tests for damage materials
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
#include "test_material_fixtures.hh"
/* -------------------------------------------------------------------------- */
#include <material_marigo.hh>
#include <material_mazars.hh>
#include <py_aka_array.hh>
#include <solid_mechanics_model.hh>
/* -------------------------------------------------------------------------- */
#include <fstream>
#include <gtest/gtest.h>
#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <type_traits>
/* -------------------------------------------------------------------------- */

using namespace akantu;

namespace py = pybind11;
using namespace py::literals;

using mat_types = ::testing::Types<
    // Traits<MaterialMarigo, 1>, Traits<MaterialMarigo, 2>,
    // Traits<MaterialMarigo, 3>,
    Traits<MaterialMazars, 1>, Traits<MaterialMazars, 2>,
    Traits<MaterialMazars, 3>>;

/*****************************************************************/

template <> void FriendMaterial<MaterialMazars<1>>::setParams() {
  K0.setDefaultValue(1e-4);
  At = 1.0;
  Bt = 5e3;
  Ac = 0.8;
  Bc = 1391.3;
  beta = 1.;
  E = 25e9;
  nu = 0.2;

  updateInternalParameters();
}

template <> void FriendMaterial<MaterialMazars<2>>::setParams() {
  K0.setDefaultValue(1e-4);
  At = 1.0;
  Bt = 5e3;
  Ac = 0.8;
  Bc = 1391.3;
  beta = 1.;
  E = 25e9;
  nu = 0.2;
  plane_stress = true;
  updateInternalParameters();
}

template <> void FriendMaterial<MaterialMazars<3>>::setParams() {
  K0.setDefaultValue(1e-4);
  At = 1.0;
  Bt = 5e3;
  Ac = 0.8;
  Bc = 1391.3;
  beta = 1.;
  E = 25e9;
  nu = 0.2;

  updateInternalParameters();
}

template <> void FriendMaterial<MaterialMazars<1>>::testComputeStress() {
  Array<Real> epsilons(1001, 1);
  Array<Real> sigmas(1001, 1);
  Array<Real> damages(1001, 1);

  for (auto && data : enumerate(epsilons)) {
    std::get<1>(data) = 2e-6 * std::get<0>(data);
  }
  Real _K0 = K0;
  py::module py_engine = py::module::import("py_mazars");

  auto kwargs_mat_params =
      py::dict("K0"_a = _K0, "At"_a = At, "Bt"_a = Bt, "Ac"_a = Ac, "Bc"_a = Bc,
               "E"_a = E, "nu"_a = nu);
  auto kwargs = py::dict("epsilons"_a = epsilons, "sigmas"_a = sigmas,
                         "damages"_a = damages);

  auto py_mazars = py_engine.attr("Mazars")(**kwargs_mat_params);
  // auto Gf_py = py_mazars.attr("compute")(**kwargs);

  Real dam = 0.;
  Real dam_ref = 0.;
  Real ehat = 0.;

  for (auto && epsilon : epsilons) {
    Matrix<Real> strain(this->spatial_dimension, this->spatial_dimension, 0.);
    Matrix<Real> sigma(this->spatial_dimension, this->spatial_dimension, 0.);
    strain(0, 0) = epsilon;

    computeStressOnQuad(strain, sigma, dam, ehat);

    Real sigma_ref;
    auto py_data =
        py_mazars.attr("compute_step")(epsilon, sigma_ref, dam_ref, false);
    std::tie(sigma_ref, dam_ref) = py::cast<std::pair<double, double>>(py_data);

    EXPECT_NEAR(sigma(0, 0), sigma_ref, 1e-5);
    EXPECT_NEAR(dam, dam_ref, 1e-10);
  }
}

template <> void FriendMaterial<MaterialMazars<2>>::testComputeStress() {
  Array<Real> epsilons(1001, 1);
  Array<Real> sigmas(1001, 1);
  Array<Real> damages(1001, 1);

  for (auto && data : enumerate(epsilons)) {
    std::get<1>(data) = 2e-6 * std::get<0>(data);
  }
  Real _K0 = K0;
  py::module py_engine = py::module::import("py_mazars");

  auto kwargs_mat_params =
      py::dict("K0"_a = _K0, "At"_a = At, "Bt"_a = Bt, "Ac"_a = Ac, "Bc"_a = Bc,
               "E"_a = E, "nu"_a = nu);
  auto kwargs = py::dict("epsilons"_a = epsilons, "sigmas"_a = sigmas,
                         "damages"_a = damages);

  auto py_mazars = py_engine.attr("Mazars")(**kwargs_mat_params);
  // auto Gf_py = py_mazars.attr("compute")(**kwargs);

  Real dam = 0.;
  Real dam_ref = 0.;
  Real ehat = 0.;

  for (auto && epsilon : epsilons) {
    Matrix<Real> strain(this->spatial_dimension, this->spatial_dimension, 0.);
    Matrix<Real> sigma(this->spatial_dimension, this->spatial_dimension, 0.);
    strain(0, 0) = epsilon;
    strain(1, 1) = -this->nu * epsilon;

    computeStressOnQuad(strain, sigma, dam, ehat);

    Real sigma_ref;
    auto py_data =
        py_mazars.attr("compute_step")(epsilon, sigma_ref, dam_ref, false);
    std::tie(sigma_ref, dam_ref) = py::cast<std::pair<double, double>>(py_data);

    EXPECT_NEAR(sigma(0, 0), sigma_ref, 1e-5);
    EXPECT_NEAR(dam, dam_ref, 1e-10);
  }
}

template <> void FriendMaterial<MaterialMazars<3>>::testComputeStress() {
  Array<Real> epsilons(1001, 1);
  Array<Real> sigmas(1001, 1);
  Array<Real> damages(1001, 1);

  for (auto && data : enumerate(epsilons)) {
    std::get<1>(data) = 2e-6 * std::get<0>(data);
  }
  Real _K0 = K0;
  py::module py_engine = py::module::import("py_mazars");

  auto kwargs_mat_params =
      py::dict("K0"_a = _K0, "At"_a = At, "Bt"_a = Bt, "Ac"_a = Ac, "Bc"_a = Bc,
               "E"_a = E, "nu"_a = nu);
  auto kwargs = py::dict("epsilons"_a = epsilons, "sigmas"_a = sigmas,
                         "damages"_a = damages);

  auto py_mazars = py_engine.attr("Mazars")(**kwargs_mat_params);
  // auto Gf_py = py_mazars.attr("compute")(**kwargs);

  Real dam = 0.;
  Real dam_ref = 0.;
  Real ehat = 0.;

  for (auto && epsilon : epsilons) {
    Matrix<Real> strain(this->spatial_dimension, this->spatial_dimension, 0.);
    Matrix<Real> sigma(this->spatial_dimension, this->spatial_dimension, 0.);
    strain(0, 0) = epsilon;
    strain(1, 1) = strain(2, 2) = -this->nu * epsilon;

    computeStressOnQuad(strain, sigma, dam, ehat);

    Real sigma_ref;
    auto py_data =
        py_mazars.attr("compute_step")(epsilon, sigma_ref, dam_ref, false);
    std::tie(sigma_ref, dam_ref) = py::cast<std::pair<double, double>>(py_data);

    EXPECT_NEAR(sigma(0, 0), sigma_ref, 1e-5);
    EXPECT_NEAR(dam, dam_ref, 1e-10);
  }
}

namespace {

template <typename T>
class TestDamageMaterialFixture : public ::TestMaterialFixture<T> {};

TYPED_TEST_SUITE(TestDamageMaterialFixture, mat_types);

TYPED_TEST(TestDamageMaterialFixture, ComputeStress) {
  this->material->testComputeStress();
}
TYPED_TEST(TestDamageMaterialFixture, DISABLED_EnergyDensity) {
  this->material->testEnergyDensity();
}
TYPED_TEST(TestDamageMaterialFixture, DISABLED_ComputeTangentModuli) {
  this->material->testComputeTangentModuli();
}
TYPED_TEST(TestDamageMaterialFixture, DISABLED_ComputeCelerity) {
  this->material->testCelerity();
}
} // namespace
/*****************************************************************/
