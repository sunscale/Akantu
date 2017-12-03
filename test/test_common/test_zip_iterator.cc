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
#include <iostream>
#include <iterator>
#include <vector>
#include <boost/range/combine.hpp>
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

template <class T> class B {
public:
  B() = default;
  B(T a) : a(a){};
  B(const B & other) = delete;
  B & operator=(const B & other) = delete;

  B(B && other)
      : a(std::move(other.a)), counter(std::move(other.counter) + 1) {}
  B & operator=(B && other) {
    a = std::move(other.a);
    counter = std::move(other.counter) + 1;
    return *this;
  }

  B & operator*=(const T & b) {
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

template <typename T> struct C {
  struct iterator {
    using reference = B<T>;
    using difference_type = void;
    using iterator_category = std::random_access_iterator_tag;
    iterator(T pos) : pos(std::move(pos)) {}

    B<T> operator*() { return B<int>(pos); }
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

template <class T>
std::ostream & operator<<(std::ostream & stream, const B<T> & b) {
  b.to_stream(stream);
  return stream;
}

// namespace std {
// template<typename T> struct iterator_traits<typename C<T>::iterator> {
//   using reference = B<T>;
//   using iterator_category = std::forward_iterator_tag;
// };
// }


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

  // ait = a.begin();
  // bit = b.begin();
  // for (; ait != aend; ++ait, ++bit) {
  //   std::cout << *ait << " " << *bit << std::endl;
  // }

  auto C1 = C<int>(0, 10);
  auto C2 = C<int>(100, 110);

  auto c1end = C1.end();
  auto c1it = C1.begin();
  auto c2it = C2.begin();

  for (;c1it != c1end; ++c1it, ++c2it) {
    std::cout << *c1it << " " << *c2it << std::endl;
  }

  for (auto && tuple : zip(C<int>(0, 10), C<int>(100, 110))) {
    std::cout << boost::get<0>(tuple) << " " << boost::get<1>(tuple) << std::endl;
  }

  return 0;
}
