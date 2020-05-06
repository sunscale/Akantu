/**
 * @file   test_facet_element_mapping.cc
 *
 * @author Dana Christen <dana.christen@gmail.com>
 *
 * @date creation: Fri May 03 2013
 * @date last modification: Tue Nov 07 2017
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
#include "aka_common.hh"
#include "aka_error.hh"
#include "mesh.hh"
#include "mesh_utils.hh"
/* -------------------------------------------------------------------------- */
#include <iostream>
#include <string>

/* -------------------------------------------------------------------------- */
using namespace akantu;
using namespace std;

int main(int argc, char * argv[]) {
  // Testing the subelement-to-element mappings
  UInt spatial_dimension(3);

  akantu::initialize(argc, argv);

  Mesh mesh(spatial_dimension, "my_mesh");
  mesh.read("./cube_physical_names.msh");

  typedef Array<std::vector<Element>> ElemToSubelemMapping;
  typedef Array<Element> SubelemToElemMapping;

  std::cout << "ELEMENT-SUBELEMENT MAPPING:" << std::endl;

  for (auto ghost_type : ghost_types) {
    std::cout << "  "
              << "Ghost type: " << ghost_type << std::endl;
    for (auto & type : mesh.elementTypes(spatial_dimension, ghost_type)) {

      const SubelemToElemMapping & subelement_to_element =
          mesh.getSubelementToElement(type, ghost_type);
      std::cout << "  "
                << "  "
                << "Element type: " << type << std::endl;

      std::cout << "  "
                << "  "
                << "  "
                << "subelement_to_element:" << std::endl;
      subelement_to_element.printself(std::cout, 8);

      for (UInt i(0); i < subelement_to_element.size(); ++i) {
        std::cout << "        ";
        for (UInt j(0); j < mesh.getNbFacetsPerElement(type); ++j) {
          if (subelement_to_element(i, j) != ElementNull) {
            std::cout << subelement_to_element(i, j);
            std::cout << ", ";
          } else {
            std::cout << "ElementNull"
                      << ", ";
          }
        }
        std::cout << "for element " << i << std::endl;
      }
    }

    for (auto & type : mesh.elementTypes(spatial_dimension - 1, ghost_type)) {
      const ElemToSubelemMapping & element_to_subelement =
          mesh.getElementToSubelement(type, ghost_type);
      std::cout << "  "
                << "  "
                << "Element type: " << type << std::endl;

      std::cout << "  "
                << "  "
                << "  "
                << "element_to_subelement:" << std::endl;
      element_to_subelement.printself(std::cout, 8);

      for (UInt i(0); i < element_to_subelement.size(); ++i) {
        const std::vector<Element> & vec = element_to_subelement(i);
        std::cout << "          ";
        std::cout << "item " << i << ": [ ";
        if (vec.size() > 0) {
          for (UInt j(0); j < vec.size(); ++j) {
            if (vec[j] != ElementNull) {
              std::cout << vec[j] << ", ";
            } else {
              std::cout << "ElementNull"
                        << ", ";
            }
          }
        } else {
          std::cout << "empty, ";
        }
        std::cout << "]"
                  << ", " << std::endl;
      }
      std::cout << std::endl;
    }
  }

  return 0;
}
