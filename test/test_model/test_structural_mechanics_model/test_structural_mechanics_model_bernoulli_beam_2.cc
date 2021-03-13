/**
 * @file   test_structural_mechanics_model_bernoulli_beam_2.cc
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
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

using namespace akantu;

/* -------------------------------------------------------------------------- */
class TestStructBernoulli2
    : public TestStructuralFixture<element_type_t<_bernoulli_beam_2>> {
  using parent = TestStructuralFixture<element_type_t<_bernoulli_beam_2>>;

public:
  void addMaterials() override {
    mat.E = 3e10;
    mat.I = 0.0025;
    mat.A = 0.01;

    this->model->addMaterial(mat);

    mat.E = 3e10;
    mat.I = 0.00128;
    mat.A = 0.01;

    this->model->addMaterial(mat);
  }

  void assignMaterials() override {
    auto & materials = this->model->getElementMaterial(parent::type);
    materials(0) = 0;
    materials(1) = 1;
  }

  void setDirichletBCs() override {
    auto boundary = this->model->getBlockedDOFs().begin(parent::ndof);
    // clang-format off
    *boundary = {true, true, true}; ++boundary;
    *boundary = {false, true, false}; ++boundary;
    *boundary = {false, true, false}; ++boundary;
    // clang-format on
  }

  void setNeumannBCs() override {
    Real M = 3600;  // Nm
    Real q = 6000; // kN/m
    Real L = 10;    // m
    auto & forces = this->model->getExternalForce();
    forces(2, 2) = -M; // moment on last node
#if 1                  // as long as integration is not available
    forces(0, 1) = -q * L / 2;
    forces(0, 2) = -q * L * L / 12;
    forces(1, 1) = -q * L / 2;
    forces(1, 2) = q * L * L / 12;
#else
    auto & group = mesh.createElementGroup("lin_force");
    group.add({type, 0, _not_ghost});
    Vector<Real> lin_force = {0, q, 0};
    // a linear force is not actually a *boundary* condition
    // it is equivalent to a volume force
    model.applyBC(BC::Neumann::FromSameDim(lin_force), group);
#endif
    forces(2, 0) = mat.E * mat.A / 18;
  }

protected:
  StructuralMaterial mat;
};

/* -------------------------------------------------------------------------- */

TEST_F(TestStructBernoulli2, TestDisplacements) {
  model->solveStep();
  auto d1 = model->getDisplacement()(1, 2);
  auto d2 = model->getDisplacement()(2, 2);
  auto d3 = model->getDisplacement()(1, 0);

  Real tol = Math::getTolerance();

  EXPECT_NEAR(d1, 5.6 / 4800, tol);  // first rotation
  EXPECT_NEAR(d2, -3.7 / 4800, tol); // second rotation
  EXPECT_NEAR(d3, 10. / 18, tol);    // axial deformation
}
