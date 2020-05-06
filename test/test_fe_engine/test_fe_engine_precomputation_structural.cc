/**
 * @file   test_fe_engine_precomputation_structural.cc
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 *
 * @date creation: Fri Jan 26 2018
 * @date last modification: Mon Feb 19 2018
 *
 * @brief  test of the structural precomputations
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
#include "fe_engine.hh"
#include "integrator_gauss.hh"
#include "shape_structural.hh"
#include "test_fe_engine_structural_fixture.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

/* -------------------------------------------------------------------------- */
// Need a special fixture for the extra normal
class TestBernoulliB3
    : public TestFEMStructuralFixture<element_type_t<_bernoulli_beam_3>> {
  using parent = TestFEMStructuralFixture<element_type_t<_bernoulli_beam_3>>;

public:
  /// Load the mesh and provide extra normal direction
  void readMesh(std::string filename) override {
    parent::readMesh(filename);
    auto & normals = this->mesh->getElementalData<Real>("extra_normal")
                         .alloc(0, dim, type, _not_ghost);
    Vector<Real> normal = {-36. / 65, -48. / 65, 5. / 13};
    normals.push_back(normal);
  }
};

/* -------------------------------------------------------------------------- */
/// Type alias
using TestBernoulliB2 =
    TestFEMStructuralFixture<element_type_t<_bernoulli_beam_2>>;
using TestDKT18 =
    TestFEMStructuralFixture<element_type_t<_discrete_kirchhoff_triangle_18>>;
/* -------------------------------------------------------------------------- */

/// Solution for 2D rotation matrices
Matrix<Real> globalToLocalRotation(Real theta) {
  auto c = std::cos(theta);
  auto s = std::sin(theta);
  return {{c, s, 0}, {-s, c, 0}, {0, 0, 1}};
}

/* -------------------------------------------------------------------------- */
TEST_F(TestBernoulliB2, PrecomputeRotations) {
  this->fem->initShapeFunctions();
  using ShapeStruct = ShapeStructural<_ek_structural>;
  auto & shape = dynamic_cast<const ShapeStruct &>(fem->getShapeFunctions());
  auto & rot = shape.getRotations(type);

  Real a = std::atan(4. / 3);
  std::vector<Real> angles = {a, -a, 0};

  Math::setTolerance(1e-15);

  for (auto && tuple : zip(make_view(rot, ndof, ndof), angles)) {
    auto rotation = std::get<0>(tuple);
    auto angle = std::get<1>(tuple);
    auto rotation_error = (rotation - globalToLocalRotation(angle)).norm<L_2>();
    EXPECT_NEAR(rotation_error, 0., Math::getTolerance());
  }
}

/* -------------------------------------------------------------------------- */
TEST_F(TestBernoulliB3, PrecomputeRotations) {
  this->fem->initShapeFunctions();
  using ShapeStruct = ShapeStructural<_ek_structural>;
  auto & shape = dynamic_cast<const ShapeStruct &>(fem->getShapeFunctions());
  auto & rot = shape.getRotations(type);

  Matrix<Real> ref = {{3. / 13, 4. / 13, 12. / 13},
                      {-4. / 5, 3. / 5, 0},
                      {-36. / 65, -48. / 65, 5. / 13}};
  Matrix<Real> solution{ndof, ndof};
  solution.block(ref, 0, 0);
  solution.block(ref, dim, dim);

  // The default tolerance is too much, really
  Math::setTolerance(1e-15);

  for (auto & rotation : make_view(rot, ndof, ndof)) {
    auto rotation_error = (rotation - solution).norm<L_2>();
    EXPECT_NEAR(rotation_error, 0., Math::getTolerance());
  }
}

/* -------------------------------------------------------------------------- */
TEST_F(TestDKT18, DISABLED_PrecomputeRotations) {
  this->fem->initShapeFunctions();
  using ShapeStruct = ShapeStructural<_ek_structural>;
  auto & shape = dynamic_cast<const ShapeStruct &>(fem->getShapeFunctions());
  auto & rot = shape.getRotations(type);

  for (auto & rotation : make_view(rot, ndof, ndof)) {
    std::cout << rotation << "\n";
  }
  std::cout.flush();
}
