/**
 * @file   test_material_cohesive_fixture.hh
 *
 * @author Nicolas Richart
 *
 * @date creation  Tue Feb 20 2018
 *
 * @brief Test the traction separations laws for cohesive elements
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
#include "solid_mechanics_model_cohesive.hh"
/* -------------------------------------------------------------------------- */
#include "test_gtest_utils.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
/* -------------------------------------------------------------------------- */

using namespace akantu;

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
  }

  void TearDown() override {
    material.reset(nullptr);
    model.reset(nullptr);
    mesh.reset(nullptr);

    openings.reset(nullptr);
    tractions.reset(nullptr);
  }

  void reset() {
    openings->resize(0);
    tractions->resize(0);
  }

  void addOpening(const Vector<Real> & direction, Real start, Real stop,
                  UInt nb_steps) {
    for (auto s : arange(nb_steps)) {
      auto opening = direction * Real(s) *(stop - start) / Real(nb_steps);
      openings->push_back(opening);
    }
    tractions->resize(openings->size());
  }

  Vector<Real> randomNormal() {
    std::random_device rd;
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis;
    Vector<Real> normal(dim);
    for (auto s : arange(dim))
      normal(s) = dis(gen);

    normal.normalize();
    return normal;
  }

  Vector<Real> applyModIOpening(Real max_opening, UInt nb_steps) {
    auto normal = randomNormal();

    this->material->insertion_stress_ = this->material->sigma_c_ * normal;

    addOpening(normal, 0., max_opening, nb_steps);

    this->material->computeTractions(*openings, normal, *tractions);

    return normal;
  }

protected:
  std::unique_ptr<Mesh> mesh;
  std::unique_ptr<SolidMechanicsModelCohesive> model;
  std::unique_ptr<Material> material;

  std::unique_ptr<Array<Real>> openings;
  std::unique_ptr<Array<Real>> tractions;
};

template <template <UInt> class Mat, UInt dim>
struct TestMaterialCohesive : public Mat<dim> {
  TestMaterialCohesive(SolidMechanicsModel & model) : Mat<dim>(model, "test") {
  }

  virtual void SetUp() {}
  virtual void resetInternal() {}

  void SetUps() {
    this->initMaterial();
    this->updateInternalParameters();

    this->SetUp();
    this->resetInternals();
  }

  void resetInternals() {
    this->resetInternal();
  }

  virtual void computeTractions(const Array<Real> & /*openings*/,
                                const Vector<Real> & /*normal*/,
                                Array<Real> & /*tractions*/) {}

  Vector<Real> insertion_stress_;
  Real sigma_c_{0};
  bool is_extrinsic{true};
};

template <template <UInt> class Mat, typename dim_>
constexpr UInt TestMaterialCohesiveFixture<Mat, dim_>::dim;
