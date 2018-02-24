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
#include "aka_iterators.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
#include <vector>
/* -------------------------------------------------------------------------- */

using namespace akantu;

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

  T a;
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
  std::vector<A<int>> a;
  std::vector<A<float>> b;
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
