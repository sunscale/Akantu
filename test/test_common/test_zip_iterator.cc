/**
 * @file   test_zip_iterator.cc
 *
 * @author Nicolas Richart
 *
 * @date creation  Fri Jul 21 2017
 *
 * @brief test the zip container and iterator
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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

class TestZipFixutre : public ::testing::Test {
  void SetUp() override {
    a.reserve(size);
    b.reserve(size);

    for (size_t i = 0; i < size; ++i) {
      a.emplace_back(i);
      b.emplace_back(i + size);
    }
  }

protected:
  size_t size{20};
  std::vector<A<int>> a;
  std::vector<A<float>> b;
};

TEST_F(TestZipFixutre, SimpleTest) {
  size_t i = 0;
  std::reference_wrapper<A<int>> a;
  std::reference_wrapper<A<float>> b;
  for (auto && pair : zip(a, b)) {
    auto && a = std::get<0>(pair);
    auto && b = std::get<1>(pair);

    EXPECT_EQ(i, a.a);
    EXPECT_EQ(0, a.copy_counter);
    EXPECT_EQ(0, a.move_counter);

    EXPECT_FLOAT_EQ(i + this->size, b.a);
    EXPECT_EQ(0, b.copy_counter);
    EXPECT_EQ(0, b.move_counter);
    ++i;
  }
}

// /* --------------------------------------------------------------------------
// */ int main() {
//   auto ait = a.begin();
//   auto bit = b.begin();
//   auto aend = a.end();
//   for (; ait != aend; ++ait, ++bit) {

//   }

//   for (auto pair : zip(a, b)) {
//     std::cout << std::get<0>(pair) << " " << std::get<1>(pair) << std::endl;
//     std::get<0>(pair) *= 10;
//   }

//   ait = a.begin();
//   bit = b.begin();
//   for (; ait != aend; ++ait, ++bit) {
//     std::cout << *ait << " " << *bit << std::endl;
//   }

//   return 0;
// }
