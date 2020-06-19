/**
 * @file   test_mesh_data.cc
 *
 * @author Dana Christen <dana.christen@gmail.com>
 *
 * @date creation: Fri May 03 2013
 * @date last modification: Wed Feb 03 2016
 *
 * @brief  Test of the MeshData class
 *
 *
 * Copyright (©) 2014-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "mesh.hh"
#include "mesh_utils.hh"
/* -------------------------------------------------------------------------- */
#include <iostream>
#include <string>

#define QUOTES(x) #x
#define ADD_QUOTES(x) QUOTES(x)

#define CAT(x, y) x##_##y
#define CONCAT(x, y) CAT(x, y)

//#define TYPE std::string
//#define VALUE1 "abc"
//#define VALUE2 "qwe"

#define ELEMENT _triangle_6
#define NAME CONCAT(TYPE, data)

/* -------------------------------------------------------------------------- */
using namespace akantu;
using namespace std;

int main() {

  std::cout << "Testing with type " << ADD_QUOTES(TYPE) << " and values "
            << ADD_QUOTES(VALUE1) << "," << ADD_QUOTES(VALUE2) << "..."
            << std::endl;
  MeshData mesh_data;

  ElementType elem_type = ELEMENT;

  const std::string name = ADD_QUOTES(NAME);

  Array<TYPE> & vec =
      mesh_data.getElementalDataArrayAlloc<TYPE>(name, elem_type);
  // XXX TO DELETE
  //  vec.copy(mesh_data.getElementalDataArrayAlloc<TYPE>(name, elem_type));

  TYPE value[2] = {VALUE1, VALUE2};

  vec.push_back(value[0]);
  vec.push_back(value[1]);

  for (UInt i(0); i < 2; i++) {
    AKANTU_DEBUG_ASSERT(vec(i) == value[i], "The Array accessed through the "
                                            "getElementDataArray method does "
                                            "not contain the right value.");
  }

  std::cout << vec << std::endl;
  std::cout << mesh_data.getTypeCode(name) << std::endl;

  return EXIT_SUCCESS;
}
