#include "material_elastic.hh"
#include "solid_mechanics_model.hh"
#include <gtest/gtest.h>
#include <random>
#include <type_traits>
/* -------------------------------------------------------------------------- */

#define TO_IMPLEMENT AKANTU_EXCEPTION("TEST TO IMPLEMENT")

using namespace akantu;

/*****************************************************************/

template <typename T> class FriendMaterial : public T {

public:
  virtual void testComputeStress() { TO_IMPLEMENT; };
  virtual void testComputeTangentModuli() { TO_IMPLEMENT; };
  virtual void testEnergyDensity() { TO_IMPLEMENT; };
  virtual void testPushWaveSpeed() { TO_IMPLEMENT; }
  virtual void testShearWaveSpeed() { TO_IMPLEMENT; }

  static inline Matrix<Real> getDeviatoricStrain(Real intensity);

  static inline Matrix<Real> getHydrostaticStrain(Real intensity);

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

/*****************************************************************/
template <typename T>
Matrix<Real> FriendMaterial<T>::getDeviatoricStrain(Real intensity) {

  Matrix<Real> dev = {{0, 1, 2}, {1, 0, 3}, {2, 3, 0}};
  dev *= intensity;
  return dev;
}

/*****************************************************************/
template <typename T>
Matrix<Real> FriendMaterial<T>::getHydrostaticStrain(Real intensity) {

  Matrix<Real> dev = {{1, 0, 0}, {0, 2, 0}, {0, 0, 3}};
  dev *= intensity;
  return dev;
}

/*****************************************************************/

template <typename T> Matrix<Real> FriendMaterial<T>::getRandomRotation3d() {
  constexpr UInt Dim = 3;
  // Seed with a real random value, if available
  std::random_device rd;
  std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis;

  Vector<Real> v1(Dim);
  Vector<Real> v2(Dim);
  Vector<Real> v3(Dim);

  for (UInt i = 0; i < Dim; ++i) {
    v1(i) = dis(gen);
    v2(i) = dis(gen);
  }

  v1.normalize();
  v2.normalize();

  auto d = v1.dot(v2);
  v2 -= v1 * d;
  v2.normalize();

  d = v1.dot(v2);
  v2 -= v1 * d;
  v2.normalize();
  
  v3.crossProduct(v2, v1);

  Matrix<Real> mat(Dim, Dim);
  for (UInt i = 0; i < Dim; ++i) {
    mat(0, i) = v1[i];
    mat(1, i) = v2[i];
    mat(2, i) = v3[i];
  }

  return mat;
}

/*****************************************************************/

template <typename T> Matrix<Real> FriendMaterial<T>::getRandomRotation2d() {

  constexpr UInt Dim = 2;
  // Seed with a real random value, if available
  std::random_device rd;
  std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis;

  Vector<Real> v1(Dim);
  Vector<Real> v2(Dim);

  for (UInt i = 0; i < Dim; ++i) {
    v1(i) = dis(gen);
    v2(i) = dis(gen);
  }

  v1.normalize();
  v2.normalize();

  auto d = v1.dot(v2);
  v2 -= v1 * d;
  v2.normalize();

  Matrix<Real> mat(Dim, Dim);
  for (UInt i = 0; i < Dim; ++i) {
    mat(0, i) = v1[i];
    mat(1, i) = v2[i];
  }

  return mat;
}

/*****************************************************************/

template <typename T> class TestMaterialFixture : public ::testing::Test {

public:
  using mat_class = typename T::mat_class;

  void SetUp() override {
    constexpr auto spatial_dimension = T::Dim;
    mesh = std::make_unique<Mesh>(spatial_dimension);
    model = std::make_unique<SolidMechanicsModel>(*mesh);
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
  std::unique_ptr<SolidMechanicsModel> model;
  std::unique_ptr<friend_class> material;
};

/*****************************************************************/
template <template <UInt> class T, UInt _Dim> struct Traits {
  using mat_class = T<_Dim>;
  static constexpr UInt Dim = _Dim;
};
/*****************************************************************/
