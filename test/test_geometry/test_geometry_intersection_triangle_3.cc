/**
 * @file   test_geometry_mesh.cc
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Fri Mar 13 2015
 * @date last modification: Fri Mar 13 2015
 *
 * @brief  Tests the interface mesh generation
 *
 * @section LICENSE
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#include "aka_common.hh"

#include "mesh_segment_intersector.hh"
#include "geom_helper_functions.hh"

#include <CGAL/Cartesian.h>

#include <iostream>

/* -------------------------------------------------------------------------- */

using namespace akantu;

typedef CGAL::Cartesian<Real> K;

/* -------------------------------------------------------------------------- */

int main (int argc, char * argv[]) {
  initialize("", argc, argv);
  debug::setDebugLevel(dblError);

  Math::setTolerance(1e-10);

  Mesh mesh(2), interface_mesh(2, "interface_mesh");
  mesh.read("test_geometry_triangle.msh");

  MeshSegmentIntersector<2, _triangle_3> intersector(mesh, interface_mesh);
  intersector.constructData();

  // Testing a segment going out of the mesh
  K::Point_3 a(0, 0.25, 0),
             b(1, 0.25, 0),
             c(0.25, 0, 0),
             d(0.25, 1, 0);

  K::Segment_3 h_interface(a, b),
               v_interface(c, d);

  std::list<K::Segment_3> interface_list;
  interface_list.push_back(h_interface);
  interface_list.push_back(v_interface);

  intersector.computeIntersectionQueryList(interface_list);

  if (interface_mesh.getNbElement(_segment_2) != 4)
    return EXIT_FAILURE;

  Vector<Real> bary(2);
  Element test;
  test.element = 0;
  test.type = _segment_2;
  
  interface_mesh.getBarycenter(test, bary);
  Real first_bary[] = {0.5, 0.25};

  if (!Math::are_vector_equal(2, bary.storage(), first_bary))
    return EXIT_FAILURE;

  // Testing a segment completely inside an element
  K::Point_3 e(0.1, 0.33, 0),
             f(0.1, 0.67, 0);
  K::Segment_3 inside_segment(e, f);
  intersector.computeIntersectionQuery(inside_segment);

  test.element = interface_mesh.getNbElement(_segment_2) - 1;
  interface_mesh.getBarycenter(test, bary);

  Real second_bary[] = {0.1, 0.5};

  if (!Math::are_vector_equal(2, bary.storage(), second_bary))
    return EXIT_FAILURE;

  finalize();
  return EXIT_SUCCESS;
}


