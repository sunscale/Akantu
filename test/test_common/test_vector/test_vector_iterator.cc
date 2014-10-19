/**
 * @file   test_vector_iterator.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Jul 30 2012
 * @date last modification: Fri Mar 21 2014
 *
 * @brief  test the iterator present in the vector class
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include <iostream>
#include <algorithm>
#include <vector>
/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_array.hh"
#include "aka_types.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

#define N 10

int main(int argc, char *argv[]) {
  typedef Array<double> RealArray;

  std::cout << "Creating a vector of matrices (2,2)" << std::endl;
  RealArray mat_vect(N, 4, 0.);

  std::cout << "Iterating on a Matrix(2,2)" << std::endl;
  RealArray::iterator<akantu::Matrix<Real> > itm;
  itm = mat_vect.begin(2, 2);
  RealArray::iterator<akantu::Matrix<Real> > endm = mat_vect.end(2, 2);

  for (; itm != endm; ++itm) {
    std::cout << *itm << std::endl;
  }

  std::cout << "Creating a vector of UInt" << std::endl;
  Array<UInt> vect(N, 1, 0.);

  std::cout << "Iterating on a UInt" << std::endl;
  Array<UInt>::iterator<UInt> it  = vect.begin();
  Array<UInt>::iterator<UInt> end = vect.end();

  std::vector<UInt> test_vect;

  for (; it != end; ++it) {
    UInt r = rand() % N;
    *it = r;
    test_vect.push_back(r);
    std::cout << *it <<std::endl;
  }

  std::cout << "Sorting" << std::endl;  
  std::sort(vect.begin(), vect.end());
  std::sort(test_vect.begin(), test_vect.end());
  
  std::vector<UInt>::iterator itv = test_vect.begin();
  for (it = vect.begin(); it != end; ++it, ++itv) {
    std::cout << *it << " " << *itv << std::endl;
    if(*it != *itv)
      AKANTU_DEBUG_ERROR("The sort of the vector gived bad results ("<< *it << " != " << *itv << ")");
  }

  return EXIT_SUCCESS;
}
