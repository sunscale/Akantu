/**
 * @file   test_mesh_iterators.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Wed Aug 16 2017
 *
 * @brief  Test the mesh iterators
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "element_group.hh"
#include "mesh.hh"
#include "mesh_iterators.hh"
#include "node_group.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char * argv[]) {
  initialize(argc, argv);

  Mesh mesh(3);
  const Mesh & cmesh = mesh;
  mesh.read("iterators_mesh.msh");

  std::cout << "ElementGroups" << std::endl;
  for (auto && element_group : mesh.iterateElementGroups()) {
    std::cout << element_group.getName() << " " << element_group.getDimension()
              << std::endl;
  }

  std::cout << "NodeGroups" << std::endl;
  for (auto && node_group : cmesh.iterateNodeGroups()) {
    std::cout << node_group.getName() << std::endl;
  }

  std::cout << "enumerate(ElementGroups)" << std::endl;
  for (auto && element_group : enumerate(mesh.iterateElementGroups())) {
    std::cout << std::get<0>(element_group) << " "
              << std::get<1>(element_group).getName() << std::endl;
  }

  // for (auto && node_group :
  //        counting(NodeGroupsIterable(cmesh))) {
  //   std::cout << std::get<0>(node_group) << " " <<
  //   std::get<1>(node_group).getName() << std::endl;
  // }

  finalize();

  return EXIT_SUCCESS;
}
