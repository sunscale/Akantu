/**
 * @file   test_array.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Nov 09 2017
 * @date last modification: Fri Jan 26 2018
 *
 * @brief  Test the arry class
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
#include "test_gtest_utils.hh"
/* -------------------------------------------------------------------------- */
#include <aka_array.hh>
#include <aka_types.hh>
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
#include <memory>
#include <typeindex>
#include <typeinfo>
/* -------------------------------------------------------------------------- */

using namespace akantu;

namespace {

class NonTrivial {
public:
  NonTrivial() = default;
  NonTrivial(int a) : a(a){};

  bool operator==(const NonTrivial & rhs) { return a == rhs.a; }
  int a{0};
};

bool operator==(const int & a, const NonTrivial & rhs) { return a == rhs.a; }

std::ostream & operator<<(std::ostream & stream, const NonTrivial & _this) {
  stream << _this.a;
  return stream;
}

/* -------------------------------------------------------------------------- */
using TestTypes = ::testing::Types<Real, UInt, NonTrivial>;
/* -------------------------------------------------------------------------- */

::testing::AssertionResult AssertType(const char * /*a_expr*/,
                                      const char * /*b_expr*/,
                                      const std::type_info & a,
                                      const std::type_info & b) {
  if (std::type_index(a) == std::type_index(b))
    return ::testing::AssertionSuccess();

  return ::testing::AssertionFailure()
         << debug::demangle(a.name()) << " != " << debug::demangle(b.name())
         << ") are different";
}

/* -------------------------------------------------------------------------- */

template <typename T> class ArrayConstructor : public ::testing::Test {
protected:
  using type = T;

  void SetUp() override { type_str = debug::demangle(typeid(T).name()); }

  template <typename... P> decltype(auto) construct(P &&... params) {
    return std::make_unique<Array<T>>(std::forward<P>(params)...);
  }

protected:
  std::string type_str;
};

TYPED_TEST_SUITE(ArrayConstructor, TestTypes);

TYPED_TEST(ArrayConstructor, ConstructDefault1) {
  auto array = this->construct();
  EXPECT_EQ(0, array->size());
  EXPECT_EQ(1, array->getNbComponent());
  EXPECT_STREQ("", array->getID().c_str());
}

TYPED_TEST(ArrayConstructor, ConstructDefault2) {
  auto array = this->construct(1000);
  EXPECT_EQ(1000, array->size());
  EXPECT_EQ(1, array->getNbComponent());
  EXPECT_STREQ("", array->getID().c_str());
}

TYPED_TEST(ArrayConstructor, ConstructDefault3) {
  auto array = this->construct(1000, 10);
  EXPECT_EQ(1000, array->size());
  EXPECT_EQ(10, array->getNbComponent());
  EXPECT_STREQ("", array->getID().c_str());
}

TYPED_TEST(ArrayConstructor, ConstructDefault4) {
  auto array = this->construct(1000, 10, "test");
  EXPECT_EQ(1000, array->size());
  EXPECT_EQ(10, array->getNbComponent());
  EXPECT_STREQ("test", array->getID().c_str());
}

TYPED_TEST(ArrayConstructor, ConstructDefault5) {
  auto array = this->construct(1000, 10, 1);
  EXPECT_EQ(1000, array->size());
  EXPECT_EQ(10, array->getNbComponent());
  EXPECT_EQ(1, array->operator()(10, 6));
  EXPECT_STREQ("", array->getID().c_str());
}

// TYPED_TEST(ArrayConstructor, ConstructDefault6) {
//   typename TestFixture::type defaultv[2] = {0, 1};

//   auto array = this->construct(1000, 2, defaultv);
//   EXPECT_EQ(1000, array->size());
//   EXPECT_EQ(2, array->getNbComponent());
//   EXPECT_EQ(1, array->operator()(10, 1));
//   EXPECT_EQ(0, array->operator()(603, 0));
//   EXPECT_STREQ("", array->getID().c_str());
// }

/* -------------------------------------------------------------------------- */
template <typename T> class ArrayFixture : public ArrayConstructor<T> {
public:
  void SetUp() override {
    ArrayConstructor<T>::SetUp();
    array = this->construct(1000, 10);
  }

  void TearDown() override { array.reset(nullptr); }

protected:
  std::unique_ptr<Array<T>> array;
};

TYPED_TEST_SUITE(ArrayFixture, TestTypes);

TYPED_TEST(ArrayFixture, Copy) {
  Array<typename TestFixture::type> copy(*this->array);

  EXPECT_EQ(1000, copy.size());
  EXPECT_EQ(10, copy.getNbComponent());
  EXPECT_NE(this->array->storage(), copy.storage());
}

TYPED_TEST(ArrayFixture, Set) {
  auto & arr = *(this->array);
  arr.set(12);
  EXPECT_EQ(12, arr(156, 5));
  EXPECT_EQ(12, arr(520, 7));
  EXPECT_EQ(12, arr(999, 9));
}

TYPED_TEST(ArrayFixture, Resize) {
  auto & arr = *(this->array);

  auto * ptr = arr.storage();

  arr.resize(0);
  EXPECT_EQ(0, arr.size());
  EXPECT_TRUE(arr.storage() == nullptr or arr.storage() == ptr);
  EXPECT_LE(0, arr.getAllocatedSize());

  arr.resize(3000);
  EXPECT_EQ(3000, arr.size());
  EXPECT_LE(3000, arr.getAllocatedSize());

  ptr = arr.storage();

  arr.resize(0);
  EXPECT_EQ(0, arr.size());
  EXPECT_TRUE(arr.storage() == nullptr or arr.storage() == ptr);
  EXPECT_LE(0, arr.getAllocatedSize());
}

TYPED_TEST(ArrayFixture, PushBack) {
  auto & arr = *(this->array);

  auto * ptr = arr.storage();

  arr.resize(0);
  EXPECT_EQ(0, arr.size());
  EXPECT_TRUE(arr.storage() == nullptr or arr.storage() == ptr);
  EXPECT_LE(0, arr.getAllocatedSize());

  arr.resize(3000);
  EXPECT_EQ(3000, arr.size());
  EXPECT_LE(3000, arr.getAllocatedSize());

  ptr = arr.storage();

  arr.resize(0);
  EXPECT_EQ(0, arr.size());
  EXPECT_TRUE(arr.storage() == nullptr or arr.storage() == ptr);
  EXPECT_LE(0, arr.getAllocatedSize());
}

TYPED_TEST(ArrayFixture, ViewVector) {
  auto && view = make_view(*this->array, 10);
  EXPECT_NO_THROW(view.begin());
  {
    auto it = view.begin();
    EXPECT_EQ(10, it->size());
    EXPECT_PRED_FORMAT2(AssertType, typeid(*it),
                        typeid(Vector<typename TestFixture::type>));
    EXPECT_PRED_FORMAT2(AssertType, typeid(it[0]),
                        typeid(VectorProxy<typename TestFixture::type>));
  }
}

TYPED_TEST(ArrayFixture, ViewMatrix) {
  {
    auto && view = make_view(*this->array, 2, 5);

    EXPECT_NO_THROW(view.begin());
    {
      auto it = view.begin();
      EXPECT_EQ(10, it->size());
      EXPECT_EQ(2, it->size(0));
      EXPECT_EQ(5, it->size(1));

      EXPECT_PRED_FORMAT2(AssertType, typeid(*it),
                          typeid(Matrix<typename TestFixture::type>));
      EXPECT_PRED_FORMAT2(AssertType, typeid(it[0]),
                          typeid(MatrixProxy<typename TestFixture::type>));
    }
  }
}

TYPED_TEST(ArrayFixture, ViewVectorWrong) {
  auto && view = make_view(*this->array, 11);
  EXPECT_THROW(view.begin(), debug::ArrayException);
}

TYPED_TEST(ArrayFixture, ViewMatrixWrong) {
  auto && view = make_view(*this->array, 3, 7);
  EXPECT_THROW(view.begin(), debug::ArrayException);
}

TYPED_TEST(ArrayFixture, ViewMatrixIter) {
  std::size_t count = 0;
  for (auto && mat : make_view(*this->array, 10, 10)) {
    EXPECT_EQ(100, mat.size());
    EXPECT_EQ(10, mat.size(0));
    EXPECT_EQ(10, mat.size(1));
    EXPECT_PRED_FORMAT2(AssertType, typeid(mat),
                        typeid(Matrix<typename TestFixture::type>));

    ++count;
  }

  EXPECT_EQ(100, count);
}

TYPED_TEST(ArrayFixture, ConstViewVector) {
  const auto & carray = *this->array;
  auto && view = make_view(carray, 10);
  EXPECT_NO_THROW(view.begin());
  {
    auto it = view.begin();
    EXPECT_EQ(10, it->size());
    EXPECT_PRED_FORMAT2(AssertType, typeid(*it),
                        typeid(Vector<typename TestFixture::type>));
    EXPECT_PRED_FORMAT2(AssertType, typeid(it[0]),
                        typeid(VectorProxy<typename TestFixture::type>));
  }
}

} // namespace
