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
#include "aka_zip.hh"
/* -------------------------------------------------------------------------- */
#include <iostream>
#include <vector>
/* -------------------------------------------------------------------------- */

using namespace akantu;

template <class T> class A {
public:
  A() = default;
  A(T a) : a(a){};
  A(const A & other) : a(other.a), counter(other.counter + 1) {}
  A & operator=(const A & other) {
    a = other.a;
    counter = other.counter + 1;
    return *this;
  }

  A & operator*=(const T & b) {
    a *= b;
    return *this;
  }
  void to_stream(std::ostream & stream) const {
    stream << a << " [" << counter << "]";
  }

private:
  T a;
  size_t counter{0};
};

template <class T>
std::ostream & operator<<(std::ostream & stream, const A<T> & a) {
  a.to_stream(stream);
  return stream;
}

/* -------------------------------------------------------------------------- */
int main() {
  std::vector<A<int>> a{1, 2, 3, 4, 5};
  const std::vector<A<float>> b{6., 7., 8., 9., 10.};

  auto aend = a.end();

  auto ait = a.begin();
  auto bit = b.begin();

  for (; ait != aend; ++ait, ++bit) {
    std::cout << *ait << " " << *bit << std::endl;
  }

  for (auto pair : zip(a, b)) {
    std::cout << std::get<0>(pair) << " " << std::get<1>(pair) << std::endl;
    std::get<0>(pair) *= 10;
  }

  ait = a.begin();
  bit = b.begin();
  for (; ait != aend; ++ait, ++bit) {
    std::cout << *ait << " " << *bit << std::endl;
  }

  return 0;
}
