/**
 * @file   test_cohesive_buildfacets_hexahedron.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date   Wed Oct 03 10:20:53 2012
 *
 * @brief  Test for cohesive elements
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
#include <iostream>
#include <limits>
#include <fstream>

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"
#include "mesh_utils.hh"
#include "solid_mechanics_model.hh"
#include "material.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char *argv[]) {
  initialize(argc, argv);

  const UInt spatial_dimension = 3;
  const ElementType type = _hexahedron_8;

  Mesh mesh(spatial_dimension);
  mesh.read("hexahedron.msh");
  Mesh mesh_facets(mesh.initMeshFacets("mesh_facets"));

  MeshUtils::buildAllFacets(mesh, mesh_facets);

  // debug::setDebugLevel(dblDump);
  // std::cout << mesh << std::endl;
  // std::cout << mesh_facets << std::endl;

  const ElementType type_facet = mesh.getFacetType(type);
  const ElementType type_subfacet = mesh.getFacetType(type_facet);
  const ElementType type_subsubfacet = mesh.getFacetType(type_subfacet);

  /* ------------------------------------------------------------------------ */
  /* Element to Subelement testing                                            */
  /* ------------------------------------------------------------------------ */

  const Array< std::vector<Element> > & el_to_subel3 = mesh_facets.getElementToSubelement(type_facet);
  const Array< std::vector<Element> > & el_to_subel2 = mesh_facets.getElementToSubelement(type_subfacet);
  const Array< std::vector<Element> > & el_to_subel1 = mesh_facets.getElementToSubelement(type_subsubfacet);

  /// build vectors for comparison
  Array<Element> hexahedron(2);
  hexahedron(0).type = type;
  hexahedron(0).element = 0;

  hexahedron(1).type = type;
  hexahedron(1).element = 3;

  Array<Element> quadrangle(4);
  quadrangle(0).type = type_facet;
  quadrangle(0).element = 1;

  quadrangle(1).type = type_facet;
  quadrangle(1).element = 2;

  quadrangle(2).type = type_facet;
  quadrangle(2).element = 7;

  quadrangle(3).type = type_facet;
  quadrangle(3).element = 11;

  Array<Element> segment(5);
  segment(0).type = type_subfacet;
  segment(0).element = 0;

  segment(1).type = type_subfacet;
  segment(1).element = 1;

  segment(2).type = type_subfacet;
  segment(2).element = 4;

  segment(3).type = type_subfacet;
  segment(3).element = 15;

  segment(4).type = type_subfacet;
  segment(4).element = 22;


  /// comparison

  for (UInt i = 0; i < hexahedron.getSize(); ++i) {
    if (hexahedron(i).type != el_to_subel3(1)[i].type ||
  	hexahedron(i).element != el_to_subel3(1)[i].element) {
      std::cout << hexahedron(i).element << " " << el_to_subel3(4)[i].element << std::endl;
      std::cout << "The two hexahedrons connected to quadrangle 1 are wrong"
  		<< std::endl;
      return EXIT_FAILURE;
    }
  }

  for (UInt i = 0; i < quadrangle.getSize(); ++i) {
    if (quadrangle(i).type != el_to_subel2(4)[i].type ||
  	quadrangle(i).element != el_to_subel2(4)[i].element) {
      std::cout << "The quadrangles connected to segment 4 are wrong"
  		<< std::endl;
      return EXIT_FAILURE;
    }
  }

  for (UInt i = 0; i < segment.getSize(); ++i) {
    if (segment(i).type != el_to_subel1(1)[i].type ||
  	segment(i).element != el_to_subel1(1)[i].element) {
      std::cout << "The segments connected to point 1 are wrong"
  		<< std::endl;
      return EXIT_FAILURE;
    }
  }


  /* ------------------------------------------------------------------------ */
  /* Subelement to Element testing                                            */
  /* ------------------------------------------------------------------------ */

  const Array<Element> & subel_to_el3 = mesh_facets.getSubelementToElement(type);
  const Array<Element> & subel_to_el2 = mesh_facets.getSubelementToElement(type_facet);
  const Array<Element> & subel_to_el1 = mesh_facets.getSubelementToElement(type_subfacet);

  /// build vectors for comparison
  Array<Element> quadrangle2(mesh.getNbFacetsPerElement(type));
  quadrangle2(0).type = type_facet;
  quadrangle2(0).element = 1;

  quadrangle2(1).type = type_facet;
  quadrangle2(1).element = 11;

  quadrangle2(2).type = type_facet;
  quadrangle2(2).element = 16;

  quadrangle2(3).type = type_facet;
  quadrangle2(3).element = 17;

  quadrangle2(4).type = type_facet;
  quadrangle2(4).element = 18;

  quadrangle2(5).type = type_facet;
  quadrangle2(5).element = 19;

  Array<Element> segment2(4);
  segment2(0).type = type_subfacet;
  segment2(0).element = 3;

  segment2(1).type = type_subfacet;
  segment2(1).element = 6;

  segment2(2).type = type_subfacet;
  segment2(2).element = 9;

  segment2(3).type = type_subfacet;
  segment2(3).element = 11;

  Array<Element> point(2);
  point(0).type = mesh.getFacetType(type_subfacet);
  point(0).element = 5;

  point(1).type = mesh.getFacetType(type_subfacet);
  point(1).element = 7;


  /// comparison

  for (UInt i = 0; i < quadrangle2.getSize(); ++i) {
    if (quadrangle2(i).type != subel_to_el3(3, i).type ||
  	quadrangle2(i).element != subel_to_el3(3, i).element) {
      std::cout << "The quadrangles connected to hexahedron 3 are wrong"
  		<< std::endl;
      return EXIT_FAILURE;
    }
  }

  for (UInt i = 0; i < segment2.getSize(); ++i) {
    if (segment2(i).type != subel_to_el2(4, i).type ||
  	segment2(i).element != subel_to_el2(4, i).element) {
      std::cout << "The segments connected to quadrangle 4 are wrong"
  		<< std::endl;
      return EXIT_FAILURE;
    }
  }

  for (UInt i = 0; i < point.getSize(); ++i) {
    if (point(i).type != subel_to_el1(11, i).type ||
  	point(i).element != subel_to_el1(11, i).element) {
      std::cout << "The points connected to segment 11 are wrong"
  		<< std::endl;
      return EXIT_FAILURE;
    }
  }


  finalize();

  std::cout << "OK: test_cohesive_buildfacets was passed!" << std::endl;

  return EXIT_SUCCESS;
}
