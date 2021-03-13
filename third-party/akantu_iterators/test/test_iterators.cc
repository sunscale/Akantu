/**
 * @file   test_zip_iterator.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jul 21 2017
 * @date last modification: Fri Dec 08 2017
 *
 * @brief  test the zip container and iterator
 *
 *
 * Copyright (©) 2016-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * akantu-iterators is free  software: you can redistribute it and/or  modify it
 * under the terms  of the  GNU Lesser  General Public  License as published by
 * the Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * akantu-iterators is  distributed in the  hope that it  will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public
 * License  for more details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with akantu-iterators. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_iterators.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
#include <vector>
/* -------------------------------------------------------------------------- */

using namespace aka;

/* -------------------------------------------------------------------------- */

// Non Trivial class that counts moves and copies
template <class T> class A {
public:
  A() = default;
  A(T a) : a(a){};
  A(const A & other)
      : a(other.a), copy_counter(other.copy_counter + 1),
        move_counter(other.move_counter) {}
  A & operator=(const A & other) {
    if (this != &other) {
      a = other.a;
      copy_counter = other.copy_counter + 1;
    }
    return *this;
  }

  A(A && other)
      : a(std::move(other.a)), copy_counter(std::move(other.copy_counter)),
        move_counter(std::move(other.move_counter) + 1) {}

  A & operator=(A && other) {
    if (this != &other) {
      a = std::move(other.a);
      copy_counter = std::move(other.copy_counter);
      move_counter = std::move(other.move_counter) + 1;
    }
    return *this;
  }

  A & operator*=(const T & b) {
    a *= b;
    return *this;
  }

  T a{};
  size_t copy_counter{0};
  size_t move_counter{0};
};

template <typename T> struct C {
  struct iterator {
    using reference = A<T>;
    using difference_type = void;
    using iterator_category = std::input_iterator_tag;
    using value_type = A<T>;
    using pointer = A<T> *;

    iterator(T pos) : pos(std::move(pos)) {}

    A<T> operator*() { return A<int>(pos); }
    bool operator!=(const iterator & other) const { return pos != other.pos; }
    bool operator==(const iterator & other) const { return pos == other.pos; }
    iterator & operator++() {
      ++pos;
      return *this;
    }
    T pos;
  };

  C(T begin_, T end_) : begin_(std::move(begin_)), end_(std::move(end_)) {}

  iterator begin() { return iterator(begin_); }
  iterator end() { return iterator(end_); }

  T begin_, end_;
};

class TestZipFixutre : public ::testing::Test {
protected:
  void SetUp() override {
    a.reserve(size);
    b.reserve(size);

    for (size_t i = 0; i < size; ++i) {
      a.emplace_back(i);
      b.emplace_back(i + size);
    }
  }

  template <typename A, typename B>
  void check(A && a, B && b, size_t pos, size_t nb_copy, size_t nb_move) {
    EXPECT_EQ(pos, a.a);
    EXPECT_EQ(nb_copy, a.copy_counter);
    EXPECT_EQ(nb_move, a.move_counter);

    EXPECT_FLOAT_EQ(pos + this->size, b.a);
    EXPECT_EQ(nb_copy, b.copy_counter);
    EXPECT_EQ(nb_move, b.move_counter);
  }

protected:
  size_t size{20};
  std::vector<A<int>> a{};
  std::vector<A<float>> b{};
};

TEST_F(TestZipFixutre, SimpleTest) {
  size_t i = 0;
  for (auto && pair : zip(this->a, this->b)) {
    this->check(std::get<0>(pair), std::get<1>(pair), i, 0, 0);
    ++i;
  }
}

TEST_F(TestZipFixutre, ConstTest) {
  size_t i = 0;
  const auto & ca = this->a;
  const auto & cb = this->b;
  for (auto && pair : zip(ca, cb)) {
    this->check(std::get<0>(pair), std::get<1>(pair), i, 0, 0);
    EXPECT_EQ(true,
              std::is_const<
                  std::remove_reference_t<decltype(std::get<0>(pair))>>::value);
    EXPECT_EQ(true,
              std::is_const<
                  std::remove_reference_t<decltype(std::get<1>(pair))>>::value);
    ++i;
  }
}

TEST_F(TestZipFixutre, MixteTest) {
  size_t i = 0;
  const auto & cb = this->b;
  for (auto && pair : zip(a, cb)) {
    this->check(std::get<0>(pair), std::get<1>(pair), i, 0, 0);
    EXPECT_EQ(false,
              std::is_const<
                  std::remove_reference_t<decltype(std::get<0>(pair))>>::value);
    EXPECT_EQ(true,
              std::is_const<
                  std::remove_reference_t<decltype(std::get<1>(pair))>>::value);
    ++i;
  }
}

TEST_F(TestZipFixutre, MoveTest) {
  size_t i = 0;
  for (auto && pair :
       zip(C<int>(0, this->size), C<int>(this->size, 2 * this->size))) {
    this->check(std::get<0>(pair), std::get<1>(pair), i, 0, 1);
    ++i;
  }
}

TEST_F(TestZipFixutre, Bidirectional) {
  auto _zip = zip(a, b);
  auto begin = _zip.begin();

  auto it = begin;
  ++it;
  EXPECT_EQ(begin, --it);

  it = begin;
  EXPECT_EQ(begin, it++);
  EXPECT_EQ(begin, --it);

  auto it2 = it = begin;
  ++it;
  ++it2;
  EXPECT_EQ(it2, it--);
  EXPECT_EQ(begin, it);
}

TEST_F(TestZipFixutre, RandomAccess) {
  auto _zip = zip(a, b);
  auto begin = _zip.begin();
  auto end = _zip.end();

  auto && val5 = begin[5];
  this->check(std::get<0>(val5), std::get<1>(val5), 5, 0, 0);

  auto && val13 = begin[13];
  this->check(std::get<0>(val13), std::get<1>(val13), 13, 0, 0);

  EXPECT_EQ(end - begin, a.size());
  auto it = ++begin;
  EXPECT_EQ(begin + 1, ++it);
  EXPECT_EQ(begin, it - 1);
}

TEST_F(TestZipFixutre, Cat) {
  size_t i = 0;
  for (auto && data : make_zip_cat(zip(a, b), zip(a, b))) {
    this->check(std::get<0>(data), std::get<1>(data), i, 0, 0);
    this->check(std::get<2>(data), std::get<3>(data), i, 0, 0);
    ++i;
  }
}

TEST(TestNamedZipFixutre, Simple) {
  std::vector<int> a{0, 10, 20, 30, 40};
  std::vector<int> b{0, 1, 2, 3, 4};

  using namespace tuple;
  for (auto && data : named_zip(get<"a"_h>() = a, get<"b"_h>() = b)) {
    auto & a = tuple::get<"a"_h>(data);
    auto & b = tuple::get<"b"_h>(data);
    b *= 10;
    EXPECT_EQ(b, a);
  }

  for (auto && data : named_zip(get<"a"_h>() = a, get<"b"_h>() = b)) {
    auto & a = tuple::get<"a"_h>(data);
    auto & b = tuple::get<"b"_h>(data);
    EXPECT_EQ(b, a);
  }
}

/* -------------------------------------------------------------------------- */
TEST(TestArangeIterator, Stop) {
  size_t ref_i = 0;
  for (auto i : arange(10)) {
    EXPECT_EQ(ref_i, i);
    ++ref_i;
  }
}

TEST(TestArangeIterator, StartStop) {
  size_t ref_i = 1;
  for (auto i : arange(1, 10)) {
    EXPECT_EQ(ref_i, i);
    ++ref_i;
  }
}

TEST(TestArangeIterator, StartStopStep) {
  size_t ref_i = 1;
  for (auto i : arange(1, 22, 2)) {
    EXPECT_EQ(ref_i, i);
    ref_i += 2;
  }
}

TEST(TestArangeIterator, StartStopStepZipped) {
  int ref_i1 = -1, ref_i2 = 1;
  for (auto && i : zip(arange(-1, -10, -1), arange(1, 18, 2))) {
    EXPECT_EQ(ref_i1, std::get<0>(i));
    EXPECT_EQ(ref_i2, std::get<1>(i));
    ref_i1 += -1;
    ref_i2 += 2;
  }
}

/* -------------------------------------------------------------------------- */
TEST(TestEnumerateIterator, SimpleTest) {
  std::vector<int> a{0, 10, 20, 30, 40};
  std::vector<int> b{0, 2, 4, 6, 8};
  for (auto && data : enumerate(a, b)) {
    EXPECT_EQ(std::get<0>(data) * 10, std::get<1>(data));
    EXPECT_EQ(std::get<0>(data) * 2, std::get<2>(data));
  }
}

/* -------------------------------------------------------------------------- */
TEST(TestTransformAdaptor, Keys) {
  std::map<std::string, int> map{
      {"1", 1}, {"2", 2}, {"3", 3}, {"3", 3}, {"4", 4}};

  char counter = '1';
  for (auto && key : make_keys_adaptor(map)) {
    EXPECT_EQ(counter, key[0]);
    ++counter;
  }
}

TEST(TestTransformAdaptor, Values) {
  std::map<std::string, int> map{
      {"1", 1}, {"2", 2}, {"3", 3}, {"3", 3}, {"4", 4}};

  int counter = 1;
  for (auto && value : make_values_adaptor(map)) {
    EXPECT_EQ(counter, value);
    ++counter;
  }
}

static int plus1(int value) { return value + 1; }

struct Plus {
  Plus(int a) : a(a) {}
  int operator()(int b) { return a + b; }

private:
  int a{0};
};

TEST(TestTransformAdaptor, Lambda) {
  auto && container = arange(10);

  for (auto && data :
       zip(container, make_transform_adaptor(container, [](auto && value) {
             return value + 1;
           }))) {
    EXPECT_EQ(std::get<0>(data) + 1, std::get<1>(data));
  }
}

TEST(TestTransformAdaptor, LambdaLambda) {
  std::map<std::string, int> map{
      {"1", 1}, {"2", 2}, {"3", 3}, {"3", 3}, {"4", 4}};

  int counter = 1;
  for (auto && data : make_transform_adaptor(
           make_values_adaptor(map), [](auto && value) { return value + 1; })) {
    EXPECT_EQ(counter + 1, data);
    ++counter;
  }

  auto && container = arange(10);

  for (auto && data :
       zip(container, make_transform_adaptor(container, [](auto && value) {
             return value + 1;
           }))) {
    EXPECT_EQ(std::get<0>(data) + 1, std::get<1>(data));
  }
}

TEST(TestTransformAdaptor, Function) {
  auto && container = arange(10);

  for (auto && data :
       zip(container, make_transform_adaptor(container, plus1))) {
    EXPECT_EQ(std::get<0>(data) + 1, std::get<1>(data));
  }
}

TEST(TestTransformAdaptor, Functor) {
  auto && container = arange(10);

  for (auto && data :
       zip(container, make_transform_adaptor(container, Plus(1)))) {
    EXPECT_EQ(std::get<0>(data) + 1, std::get<1>(data));
  }
}

/* -------------------------------------------------------------------------- */
TEST(TestFilteredIterator, Simple) {
  std::vector<int> values{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  std::vector<int> filter_{1, 3, 4, 10, 8, 6};
  for (auto && data : zip(filter_, filter(filter_, values))) {
    EXPECT_EQ(std::get<0>(data), std::get<1>(data));
  }
}

/* -------------------------------------------------------------------------- */
TEST(TestFilteredIterator, Temporary) {
  std::vector<int> filter_{1, 3, 4, 10, 8, 6};
  for (auto && data :
       zip(filter_, filter(filter_, std::vector<int>{0, 1, 2, 3, 4, 5, 6, 7, 8,
                                                     9, 10}))) {
    EXPECT_EQ(std::get<0>(data), std::get<1>(data));
  }
}

/* -------------------------------------------------------------------------- */
TEST(TestConcatenateIterator, SimpleTest) {
  for (auto && data : zip(arange(0, 13), concat(arange(0, 5), arange(5, 10),
                                                arange(10, 13)))) {
    EXPECT_EQ(std::get<0>(data), std::get<1>(data));
  }
}
