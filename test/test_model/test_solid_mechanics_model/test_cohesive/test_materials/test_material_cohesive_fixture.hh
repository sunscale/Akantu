/**
 * @file   test_material_cohesive_fixture.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Feb 21 2018
 *
 * @brief  Test the traction separations laws for cohesive elements
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
#include "solid_mechanics_model_cohesive.hh"
/* -------------------------------------------------------------------------- */
#include "test_gtest_utils.hh"
/* -------------------------------------------------------------------------- */
#include <fstream>
#include <gtest/gtest.h>
/* -------------------------------------------------------------------------- */

using namespace akantu;

//#define debug_

/* -------------------------------------------------------------------------- */
template <template <UInt> class Mat, typename dim_>
class TestMaterialCohesiveFixture : public ::testing::Test {
public:
  static constexpr UInt dim = dim_::value;
  using Material = Mat<dim>;

  void SetUp() override {
    mesh = std::make_unique<Mesh>(dim);
    model = std::make_unique<SolidMechanicsModelCohesive>(*mesh);
    material = std::make_unique<Material>(*model);

    material->SetUps();
    openings = std::make_unique<Array<Real>>(0, dim);
    tractions = std::make_unique<Array<Real>>(0, dim);
    reset();

    gen.seed(::testing::GTEST_FLAG(random_seed));

    normal = getRandomNormal();
    tangents = getRandomTangents();
  }

  void TearDown() override {
    material.reset(nullptr);
    model.reset(nullptr);
    mesh.reset(nullptr);

    openings.reset(nullptr);
    tractions.reset(nullptr);
  }

  void reset() {
    openings->resize(1);
    tractions->resize(1);

    openings->clear();
    tractions->clear();
  }

  /* ------------------------------------------------------------------------ */
  void addOpening(const Vector<Real> & direction, Real start, Real stop,
                  UInt nb_steps) {
    for (auto s : arange(nb_steps)) {
      auto opening =
          direction * (start + (stop - start) / Real(nb_steps) * Real(s + 1));
      openings->push_back(opening);
    }
    tractions->resize(openings->size());
  }

  /* ------------------------------------------------------------------------ */
  Vector<Real> getRandomVector() {
    std::uniform_real_distribution<Real> dis;
    Vector<Real> vector(dim);
    for (auto s : arange(dim))
      vector(s) = dis(gen);
    return vector;
  }

  Vector<Real> getRandomNormal() {
    auto normal = getRandomVector();
    normal.normalize();
#if defined(debug_)
    normal.set(0.);
    normal(0) = 1.;
#endif
    return normal;
  }

  Matrix<Real> getRandomTangents() {
    auto dim = normal.size();
    Matrix<Real> tangent(dim, dim - 1);
    if (dim == 2) {
      Math::normal2(normal.storage(), tangent(0).storage());
    }

    if (dim == 3) {
      auto v = getRandomVector();
      tangent(0) = (v - v.dot(normal) * normal).normalize();
      Math::normal3(normal.storage(), tangent(0).storage(),
                    tangent(1).storage());
    }

#if defined(debug_)
    if (dim == 2)
      tangent(0) = Vector<Real>{0., 1};
    if (dim == 3)
      tangent = Matrix<Real>{{0., 0.}, {1., 0.}, {0., 1.}};
#endif

    return tangent;
  }

  /* ------------------------------------------------------------------------ */
  void output_csv() {
    const ::testing::TestInfo * const test_info =
        ::testing::UnitTest::GetInstance()->current_test_info();

    std::ofstream cout(std::string(test_info->name()) + ".csv");
    auto print_vect_name = [&](auto name) {
      for (auto s : arange(dim)) {
        if (s != 0) {
          cout << ", ";
        }
        cout << name << "_" << s;
      }
    };
    auto print_vect = [&](const auto & vect) {
      cout << vect.dot(normal);
      if (dim > 1)
        cout << ", " << vect.dot(tangents(0));
      if (dim > 2)
        cout << ", " << vect.dot(tangents(1));
    };

    cout << "delta, ";
    print_vect_name("opening");
    cout << ", ";
    print_vect_name("traction");
    cout << std::endl;

    for (auto && data : zip(make_view(*this->openings, this->dim),
                            make_view(*this->tractions, this->dim))) {
      const auto & opening = std::get<0>(data);
      auto & traction = std::get<1>(data);

      cout << this->material->delta(opening, normal) << ", ";
      print_vect(opening);
      cout << ", ";
      print_vect(traction);
      cout << std::endl;
    }
  }

  /* ------------------------------------------------------------------------ */
  Real dissipated() {
    Vector<Real> prev_opening(dim, 0.);
    Vector<Real> prev_traction(dim, 0.);

    Real etot = 0.;
    Real erev = 0.;
    for (auto && data : zip(make_view(*this->openings, this->dim),
                            make_view(*this->tractions, this->dim))) {
      const auto & opening = std::get<0>(data);
      const auto & traction = std::get<1>(data);

      etot += (opening - prev_opening).dot(traction + prev_traction) / 2.;
      erev = traction.dot(opening) / 2.;

      prev_opening = opening;
      prev_traction = traction;
    }

    return etot - erev;
  }

  /* ------------------------------------------------------------------------ */
  void checkModeI(Real max_opening, Real expected_dissipated) {
    this->material->insertion_stress_ = this->material->sigma_c_ * normal;
    addOpening(normal, 0., max_opening, 100);

    this->material->computeTractions(*openings, normal, *tractions);

    for (auto && data : zip(make_view(*this->openings, this->dim),
                            make_view(*this->tractions, this->dim))) {
      const auto & opening = std::get<0>(data);
      auto & traction = std::get<1>(data);
      auto T = traction.dot(normal);
      EXPECT_NEAR(0, (traction - T * normal).norm(), 1e-9);

      auto T_expected =
          this->material->tractionModeI(opening, normal).dot(normal);
      EXPECT_NEAR(T_expected, T, 1e-9);
    }

    EXPECT_NEAR(expected_dissipated, dissipated(), 1e-5);
    this->output_csv();
  }

  /* ------------------------------------------------------------------------ */
  void checkModeII(Real max_opening) {
    if (this->dim == 1) {
      SUCCEED();
      return;
    }

    std::uniform_real_distribution<Real> dis;

    auto direction = Vector<Real>(tangents(0));
    auto alpha = dis(gen) + 0.1;
    auto beta = dis(gen) + 0.2;
#ifndef debug_
    direction = alpha * Vector<Real>(tangents(0));
    if (dim > 2)
      direction += beta * Vector<Real>(tangents(1));

    direction = direction.normalize();
#endif

    beta = this->material->get("beta");
    this->material->insertion_stress_ =
        beta * this->material->sigma_c_ * direction;

    addOpening(direction, 0., max_opening, 100);

    this->material->computeTractions(*openings, normal, *tractions);

    for (auto && data : zip(make_view(*this->openings, this->dim),
                            make_view(*this->tractions, this->dim))) {
      const auto & opening = std::get<0>(data);
      const auto & traction = std::get<1>(data);

      // In ModeII normal traction should be 0
      ASSERT_NEAR(0, traction.dot(normal), 1e-9);
      // Normal opening is null
      ASSERT_NEAR(0, opening.dot(normal), 1e-16);

      auto T = traction.dot(direction);
      auto T_expected =
          this->material->tractionModeII(opening, normal).dot(direction);

      EXPECT_NEAR(T_expected, T, 1e-9);
    }

    // EXPECT_NEAR(expected_dissipated, dissipated(), 1e-5);
    this->output_csv();
  }

protected:
  Vector<Real> normal;
  Matrix<Real> tangents;
  std::unique_ptr<Mesh> mesh;
  std::unique_ptr<SolidMechanicsModelCohesive> model;
  std::unique_ptr<Material> material;

  std::unique_ptr<Array<Real>> openings;
  std::unique_ptr<Array<Real>> tractions;

  std::mt19937 gen;
};

template <template <UInt> class Mat, UInt dim>
struct TestMaterialCohesive : public Mat<dim> {
  TestMaterialCohesive(SolidMechanicsModel & model)
      : Mat<dim>(model, "test"), insertion_stress_(dim, 0.) {}

  virtual void SetUp() {}
  virtual void resetInternal() {}

  void SetUps() {
    this->initMaterial();
    this->SetUp();
    this->updateInternalParameters();
    this->resetInternals();
  }

  void resetInternals() { this->resetInternal(); }

  virtual void computeTractions(Array<Real> & /*openings*/,
                                const Vector<Real> & /*normal*/,
                                Array<Real> & /*tractions*/) {}

  Vector<Real> insertion_stress_;
  Real sigma_c_{0};
  bool is_extrinsic{true};
};

template <template <UInt> class Mat, typename dim_>
constexpr UInt TestMaterialCohesiveFixture<Mat, dim_>::dim;
