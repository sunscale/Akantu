/**
 * @file   test_structural_mechanics_model_bernoulli_beam_3.cc
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 *
 * @date creation: Fri Jul 15 2011
 * @date last modification: Fri Feb 09 2018
 *
 * @brief  Computation of the analytical exemple 1.1 in the TGC vol 6
 *
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "test_structural_mechanics_model_fixture.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
/* -------------------------------------------------------------------------- */

using namespace akantu;

class TestStructBernoulli3
    : public TestStructuralFixture<element_type_t<_bernoulli_beam_3>> {
  using parent = TestStructuralFixture<element_type_t<_bernoulli_beam_3>>;

public:
  void readMesh(std::string filename) override {
    parent::readMesh(filename);
    auto & normals =
        this->mesh->getElementalData<Real>("extra_normal")
            .alloc(0, parent::spatial_dimension, parent::type, _not_ghost);
    Vector<Real> normal = {0, 0, 1};
    normals.push_back(normal);
    normal = {0, 0, 1};
    normals.push_back(normal);
  }

  void addMaterials() override {
    StructuralMaterial mat;
    mat.E = 1;
    mat.Iz = 1;
    mat.Iy = 1;
    mat.A = 1;
    mat.GJ = 1;
    this->model->addMaterial(mat);
  }

  void setDirichlets() override {
    // Boundary conditions (blocking all DOFs of nodes 2 & 3)
    auto boundary = ++this->model->getBlockedDOFs().begin(parent::ndof);
    // clang-format off
    *boundary = {true, true, true, true, true, true}; ++boundary;
    *boundary = {true, true, true, true, true, true}; ++boundary;
    // clang-format on
  }

  void setNeumanns() override {
    // Forces
    Real P = 1; // N
    auto & forces = this->model->getExternalForce();
    forces.zero();
    forces(0, 2) = -P; // vertical force on first node
  }

  void assignMaterials() override {
    model->getElementMaterial(parent::type).set(0);
  }
};

/* -------------------------------------------------------------------------- */

TEST_F(TestStructBernoulli3, TestDisplacements) {
  model->solveStep();
  auto vz = model->getDisplacement()(0, 2);
  auto thy = model->getDisplacement()(0, 4);
  auto thx = model->getDisplacement()(0, 3);
  auto thz = model->getDisplacement()(0, 5);

  Real tol = Math::getTolerance();

  EXPECT_NEAR(vz, -5. / 48, tol);
  EXPECT_NEAR(thy, -std::sqrt(2) / 8, tol);
  EXPECT_NEAR(thz, 0, tol);
  EXPECT_NEAR(thx, 0, tol);
}
