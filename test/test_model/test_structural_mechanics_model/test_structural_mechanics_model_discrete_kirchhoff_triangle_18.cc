/**
 * @file   test_structural_mechanics_model_discrete_kirchhoff_triangle_18.cc
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 *
 * @date creation: Fri Jul 15 2011
 * @date last modification: Wed Feb 21 2018
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
#include "sparse_matrix.hh"
#include "test_structural_mechanics_model_fixture.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>

using namespace akantu;

/* -------------------------------------------------------------------------- */
class TestStructDKT18 : public TestStructuralFixture<
                            element_type_t<_discrete_kirchhoff_triangle_18>> {
  using parent =
      TestStructuralFixture<element_type_t<_discrete_kirchhoff_triangle_18>>;

public:
  void addMaterials() override {
    mat.E = 1;
    mat.t = 1;
    mat.nu = 0.3;

    this->model->addMaterial(mat);
  }

  void assignMaterials() override {
    auto & materials = this->model->getElementMaterial(parent::type);
    std::fill(materials.begin(), materials.end(), 0);
  }

  void setDirichletBCs() override {
    this->model->getBlockedDOFs().set(true);
    auto center_node = this->model->getBlockedDOFs().end(parent::ndof) - 1;
    *center_node = {false, false, false, false, false, true};

    this->model->getDisplacement().zero();
    auto disp = ++this->model->getDisplacement().begin(parent::ndof);

    // Displacement field from Batoz Vol. 2 p. 392
    // with theta to beta correction from discrete Kirchhoff condition
    // see also the master thesis of Michael Lozano

    // clang-format off

    // This displacement field tests membrane and bending modes
    *disp = {40, 20, -800 , -20, 40, 0}; ++disp;
    *disp = {50, 40, -1400, -40, 50, 0}; ++disp;
    *disp = {10, 20, -200 , -20, 10, 0}; ++disp;

    // This displacement tests the bending mode
    // *disp = {0, 0, -800 , -20, 40, 0}; ++disp;
    // *disp = {0, 0, -1400, -40, 50, 0}; ++disp;
    // *disp = {0, 0, -200 , -20, 10, 0}; ++disp;

    // This displacement tests the membrane mode
    // *disp = {40, 20, 0, 0, 0, 0}; ++disp;
    // *disp = {50, 40, 0, 0, 0, 0}; ++disp;
    // *disp = {10, 20, 0, 0, 0, 0}; ++disp;

    // clang-format on
  }

  void setNeumannBCs() override {}

protected:
  StructuralMaterial mat;
};

/* -------------------------------------------------------------------------- */

// Batoz Vol 2. patch test, ISBN 2-86601-259-3
TEST_F(TestStructDKT18, TestDisplacements) {
  model->solveStep();
  Vector<Real> solution = {22.5, 22.5, -337.5, -22.5, 22.5, 0};
  auto nb_nodes = this->model->getDisplacement().size();

  Vector<Real> center_node_disp =
      model->getDisplacement().begin(solution.size())[nb_nodes - 1];

  auto error = solution - center_node_disp;

  EXPECT_NEAR(error.norm<L_2>(), 0., 1e-12);
}
