/**
 * @file   test_material_fixtures.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date creation: Fri Nov 17 2017
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Fixture for material tests
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
#include "material_elastic.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
#include <random>
#include <type_traits>
/* -------------------------------------------------------------------------- */

using namespace akantu;

/* -------------------------------------------------------------------------- */
template <typename T> class FriendMaterial : public T {
public:
  ~FriendMaterial() = default;

  virtual void testComputeStress() { AKANTU_TO_IMPLEMENT(); };
  virtual void testComputeTangentModuli() { AKANTU_TO_IMPLEMENT(); };
  virtual void testEnergyDensity() { AKANTU_TO_IMPLEMENT(); };
  virtual void testCelerity() { AKANTU_TO_IMPLEMENT(); }

  static inline Matrix<Real> getDeviatoricStrain(Real intensity);
  static inline Matrix<Real> getHydrostaticStrain(Real intensity);
  static inline Matrix<Real> getComposedStrain(Real intensity);

  static inline const Matrix<Real>
  reverseRotation(Matrix<Real> mat, Matrix<Real> rotation_matrix) {
    Matrix<Real> R = rotation_matrix;
    Matrix<Real> m2 = mat * R;
    Matrix<Real> m1 = R.transpose();
    return m1 * m2;
  };

  static inline const Matrix<Real>
  applyRotation(Matrix<Real> mat, const Matrix<Real> rotation_matrix) {
    Matrix<Real> R = rotation_matrix;
    Matrix<Real> m2 = mat * R.transpose();
    Matrix<Real> m1 = R;
    return m1 * m2;
  };

  static inline Matrix<Real> getRandomRotation3d();
  static inline Matrix<Real> getRandomRotation2d();

  using T::T;
};

/* -------------------------------------------------------------------------- */
template <typename T>
Matrix<Real> FriendMaterial<T>::getDeviatoricStrain(Real intensity) {
  Matrix<Real> dev = {{0., 1., 2.}, {1., 0., 3.}, {2., 3., 0.}};
  dev *= intensity;
  return dev;
}

/* -------------------------------------------------------------------------- */
template <typename T>
Matrix<Real> FriendMaterial<T>::getHydrostaticStrain(Real intensity) {
  Matrix<Real> dev = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  dev *= intensity;
  return dev;
}

/* -------------------------------------------------------------------------- */

template <typename T>
Matrix<Real> FriendMaterial<T>::getComposedStrain(Real intensity) {
  Matrix<Real> s = FriendMaterial<T>::getHydrostaticStrain(intensity) +
                   FriendMaterial<T>::getDeviatoricStrain(intensity);
  s *= intensity;
  return s;
}

/* -------------------------------------------------------------------------- */
template <typename T> Matrix<Real> FriendMaterial<T>::getRandomRotation3d() {
  constexpr UInt Dim = 3;
  // Seed with a real random value, if available
  std::random_device rd;
  std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<Real> dis;

  Vector<Real> v1(Dim);
  Vector<Real> v2(Dim);
  Vector<Real> v3(Dim);

  for (UInt i = 0; i < Dim; ++i) {
    v1(i) = dis(gen);
    v2(i) = dis(gen);
  }

  v1.normalize();

  auto d = v1.dot(v2);
  v2 -= v1 * d;
  v2.normalize();

  v3.crossProduct(v1, v2);

  Matrix<Real> mat(Dim, Dim);
  for (UInt i = 0; i < Dim; ++i) {
    mat(0, i) = v1[i];
    mat(1, i) = v2[i];
    mat(2, i) = v3[i];
  }

  return mat;
}

/* -------------------------------------------------------------------------- */
template <typename T> Matrix<Real> FriendMaterial<T>::getRandomRotation2d() {
  constexpr UInt Dim = 2;
  // Seed with a real random value, if available
  std::random_device rd;
  std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis;

  // v1 and v2 are such that they form a right-hand basis with prescribed v3,
  // it's need (at least) for 2d linear elastic anisotropic and (orthotropic)
  // materials to rotate the tangent stiffness

  Vector<Real> v1(3);
  Vector<Real> v2(3);
  Vector<Real> v3 = {0, 0, 1};

  for (UInt i = 0; i < Dim; ++i) {
    v1(i) = dis(gen);
    // v2(i) = dis(gen);
  }

  v1.normalize();
  v2.crossProduct(v3, v1);

  Matrix<Real> mat(Dim, Dim);
  for (UInt i = 0; i < Dim; ++i) {
    mat(0, i) = v1[i];
    mat(1, i) = v2[i];
  }

  return mat;
}

/* -------------------------------------------------------------------------- */
template <typename T, class Model>
class TestMaterialBaseFixture : public ::testing::Test {
public:
  using mat_class = typename T::mat_class;

  void SetUp() override {
    constexpr auto spatial_dimension = T::Dim;
    mesh = std::make_unique<Mesh>(spatial_dimension);
    model = std::make_unique<Model>(*mesh);
    material = std::make_unique<friend_class>(*model);
  }

  void TearDown() override {
    material.reset(nullptr);
    model.reset(nullptr);
    mesh.reset(nullptr);
  }

  using friend_class = FriendMaterial<mat_class>;

protected:
  std::unique_ptr<Mesh> mesh;
  std::unique_ptr<Model> model;
  std::unique_ptr<friend_class> material;
};

/* -------------------------------------------------------------------------- */
template <template <UInt> class T, UInt _Dim> struct Traits {
  using mat_class = T<_Dim>;
  static constexpr UInt Dim = _Dim;
};

/* -------------------------------------------------------------------------- */
template <typename T>
using TestMaterialFixture = TestMaterialBaseFixture<T, SolidMechanicsModel>;

/* -------------------------------------------------------------------------- */
