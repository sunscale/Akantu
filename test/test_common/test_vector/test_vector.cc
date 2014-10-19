/**
 * @file   test_vector.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Sep 03 2010
 * @date last modification: Thu Mar 27 2014
 *
 * @brief  test of the vector class
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
#include <cstdlib>
#include <iostream>
/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_types.hh"
#include "aka_array.hh"

/* -------------------------------------------------------------------------- */

int main(int argc, char *argv[]) {
  int def_value[3];
  def_value[0] = 10;
  def_value[1] = 20;
  def_value[2] = 30;

  std::cout << "Creating a vector" << std::endl;
  akantu::Array<int> int_vect(10, 3, def_value, "test1");

  for (unsigned int i = 5; i < int_vect.getSize(); ++i) {
    for (unsigned int j = 0; j < int_vect.getNbComponent(); ++j) {
      int_vect.storage()[i*int_vect.getNbComponent() + j] = def_value[j]*10;
    }
  }

  std::cerr << int_vect;

  int new_elem[3];
  new_elem[0] = 1;
  new_elem[1] = 2;
  new_elem[2] = 3;
  std::cout << "Testing push_back" << std::endl;
  int_vect.push_back(new_elem);

  int_vect.push_back(200);

  int_vect.erase(0);

  std::cerr << int_vect;
  akantu::Array<int> int_vect0(0, 3, def_value, "test2");
  std::cerr << int_vect0;
  int_vect0.push_back(new_elem);
  std::cerr << int_vect0;

  return EXIT_SUCCESS;
}
