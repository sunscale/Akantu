/**
 * @file   test_fe_engine_precomputation_structural.cc
 *
 * @author Lucas Frérot
 *
 * @date creation: Fri Feb 26 2018
 *
 * @brief  test of the structural precomputations
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include "fe_engine.hh"
#include "integrator_gauss.hh"
#include "shape_structural.hh"
#include "test_fe_engine_structural_fixture.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/// Testing if two collections of floats are equal
template <typename T, typename U>
::testing::AssertionResult array_equal(const T & a, const U & b) {
  if (a.size() != b.size())
    return ::testing::AssertionFailure() << "vector sizes are not equal";
  auto size = a.size();
  if (Math::are_vector_equal(size, a.storage(), b.storage()))
    return ::testing::AssertionSuccess() << "vectors are equal";
  else
    return ::testing::AssertionFailure() << "vectors are not equal";
}

/* -------------------------------------------------------------------------- */
// Need a special fixture for the extra normal
class TestBernoulliB3
    : public TestFEMStructuralFixture<element_type_t<_bernoulli_beam_3>> {
  using parent = TestFEMStructuralFixture<element_type_t<_bernoulli_beam_3>>;
public:

  /// Load the mesh and provide extra normal direction
  void readMesh(std::string filename) override {
    parent::readMesh(filename);
    auto & normals = this->mesh->registerData<Real>("extra_normal")
                         .alloc(0, dim, type, _not_ghost);
    Vector<Real> normal = {-36. / 65, -48. / 65, 5. / 13};
    normals.push_back(normal);
  }
};

/* -------------------------------------------------------------------------- */
/// Type alias
using TestBernoulliB2 =
    TestFEMStructuralFixture<element_type_t<_bernoulli_beam_2>>;
/* -------------------------------------------------------------------------- */

/// Solution for 2D rotation matrices
Matrix<Real> globalToLocalRotation(Real theta) {
  auto c = std::cos(theta);
  auto s = std::sin(theta);
  return {{c, s, 0}, {-s, c, 0}, {0, 0, 1}};
}

/* -------------------------------------------------------------------------- */
TEST_F(TestBernoulliB2, PrecomputeRotations) {
  using ShapeStruct = ShapeStructural<_ek_structural>;
  auto & shape = dynamic_cast<const ShapeStruct &>(fem->getShapeFunctions());
  auto & rot = shape.getRotations(type);

  Real a = std::atan(4. / 3);
  std::vector<Real> angles = {a, -a, 0};

  for (auto && tuple : zip(make_view(rot, ndof, ndof), angles)) {
    auto rotation = std::get<0>(tuple);
    auto angle = std::get<1>(tuple);

    EXPECT_TRUE(array_equal(rotation, globalToLocalRotation(angle)));
  }
}

/* -------------------------------------------------------------------------- */
TEST_F(TestBernoulliB3, PrecomputeRotations) {
  using ShapeStruct = ShapeStructural<_ek_structural>;
  auto & shape = dynamic_cast<const ShapeStruct &>(fem->getShapeFunctions());
  auto & rot = shape.getRotations(type);

  Matrix<Real> ref = {{3. / 13, 4. / 13, 12. / 13},
                      {-4. / 5, 3. / 5, 0},
                      {-36. / 65, -48. / 65, 5. / 13}};
  Matrix<Real> solution{ndof, ndof};
  solution.block(ref, 0, 0);
  solution.block(ref, dim, dim);

  for (auto & rotation : make_view(rot, ndof, ndof)) {
    EXPECT_TRUE(array_equal(rotation, solution));
  }
}
