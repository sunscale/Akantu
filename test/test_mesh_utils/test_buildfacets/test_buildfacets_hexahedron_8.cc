/**
 * @file   test_buildfacets_hexahedron_8.cc
 *
 * @author Mauro Corrado <mauro.corrado@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Tue May 08 2012
 * @date last modification: Thu Nov 09 2017
 *
 * @brief  Test to check the building of the facets. Mesh with hexahedrons
 *
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
#include <fstream>
#include <iostream>
#include <limits>

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_utils.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char * argv[]) {
  initialize(argc, argv);

  const UInt spatial_dimension = 3;
  const ElementType type = _hexahedron_8;

  Mesh mesh(spatial_dimension);
  mesh.read("hexahedron_8.msh");
  Mesh & mesh_facets = mesh.initMeshFacets("mesh_facets");

  // debug::setDebugLevel(dblDump);
  // std::cout << mesh << std::endl;
  // std::cout << mesh_facets << std::endl;

  const ElementType type_facet = mesh.getFacetType(type);
  const ElementType type_subfacet = mesh.getFacetType(type_facet);
  const ElementType type_subsubfacet = mesh.getFacetType(type_subfacet);

  /* ------------------------------------------------------------------------ */
  /* Element to Subelement testing                                            */
  /* ------------------------------------------------------------------------ */

  const Array<std::vector<Element>> & el_to_subel3 =
      mesh_facets.getElementToSubelement(type_facet);
  const Array<std::vector<Element>> & el_to_subel2 =
      mesh_facets.getElementToSubelement(type_subfacet);
  const Array<std::vector<Element>> & el_to_subel1 =
      mesh_facets.getElementToSubelement(type_subsubfacet);

  std::cout << "ElementToSubelement3" << std::endl;
  for (UInt i = 0; i < el_to_subel3.size(); ++i) {
    std::cout << type_facet << " " << i << " connected to ";
    for (UInt j = 0; j < 2; ++j) {
      std::cout << el_to_subel3(i)[j].type << " " << el_to_subel3(i)[j].element
                << ", ";
    }
    std::cout << " " << std::endl;
  }

  std::cout << "ElementToSubelement2" << std::endl;
  for (UInt i = 0; i < el_to_subel2.size(); ++i) {
    std::cout << type_subfacet << " " << i << " connected to ";
    for (UInt j = 0; j < el_to_subel2(i).size(); ++j) {
      std::cout << el_to_subel2(i)[j].type << " " << el_to_subel2(i)[j].element
                << ", ";
    }
    std::cout << " " << std::endl;
  }

  std::cout << "ElementToSubelement1" << std::endl;
  for (UInt i = 0; i < el_to_subel1.size(); ++i) {
    std::cout << type_subsubfacet << " " << i << " connected to ";
    for (UInt j = 0; j < el_to_subel1(i).size(); ++j) {
      std::cout << el_to_subel1(i)[j].type << " " << el_to_subel1(i)[j].element
                << ", ";
    }
    std::cout << " " << std::endl;
  }

  /* ------------------------------------------------------------------------ */
  /* Subelement to Element testing                                            */
  /* ------------------------------------------------------------------------ */

  const Array<Element> & subel_to_el3 =
      mesh_facets.getSubelementToElement(type);
  const Array<Element> & subel_to_el2 =
      mesh_facets.getSubelementToElement(type_facet);
  const Array<Element> & subel_to_el1 =
      mesh_facets.getSubelementToElement(type_subfacet);

  std::cout << " " << std::endl;
  std::cout << "SubelementToElement3" << std::endl;
  for (UInt i = 0; i < subel_to_el3.size(); ++i) {
    std::cout << type << " " << i << " connected to ";
    for (UInt j = 0; j < 6; ++j) {
      std::cout << subel_to_el3(i, j).type << " " << subel_to_el3(i, j).element
                << ", ";
    }
    std::cout << " " << std::endl;
  }

  std::cout << "SubelementToElement2" << std::endl;
  for (UInt i = 0; i < subel_to_el2.size(); ++i) {
    std::cout << type_facet << " " << i << " connected to ";
    for (UInt j = 0; j < 4; ++j) {
      std::cout << subel_to_el2(i, j).type << " " << subel_to_el2(i, j).element
                << ", ";
    }
    std::cout << " " << std::endl;
  }

  std::cout << "SubelementToElement1" << std::endl;
  for (UInt i = 0; i < subel_to_el1.size(); ++i) {
    std::cout << type_subfacet << " " << i << " connected to ";
    for (UInt j = 0; j < 2; ++j) {
      std::cout << subel_to_el1(i, j).type << " " << subel_to_el1(i, j).element
                << ", ";
    }
    std::cout << " " << std::endl;
  }

  finalize();

  return EXIT_SUCCESS;
}
