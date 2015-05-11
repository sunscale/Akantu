/**
 * @file   test_geometry_intersection_tetrahedron_4.cc
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Thu Mar 26 2015
 * @date last modification: Thu Mar 26 2015
 *
 * @brief  Tests the intersection module with _tetrahedron_4 elements
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

#include <CGAL/Cartesian.h>

#include <iostream>

/* -------------------------------------------------------------------------- */

using namespace akantu;

typedef CGAL::Cartesian<Real> K;
typedef K::Point_3 Point;
typedef K::Segment_3 Segment;

/* -------------------------------------------------------------------------- */

int main (int argc, char * argv[]) {
  initialize("", argc, argv);
  debug::setDebugLevel(dblError);

  Mesh mesh(3), interface_mesh(3, "interface_mesh");
  mesh.read("test_geometry_tetrahedron.msh");

  MeshSegmentIntersector<3, _tetrahedron_4> intersector(mesh, interface_mesh);
  intersector.constructData();

  // Testing a segment going through the cube
  Point point(1., 1., 1.);
  Segment segment(CGAL::ORIGIN, point);

  intersector.computeIntersectionQuery(segment);

  if (interface_mesh.getNbElement(_segment_2) != 2)
    return EXIT_FAILURE;

  Element test;
  test.type = _segment_2;
  test.element = 0;

  Vector<Real> bary(2);
  interface_mesh.getBarycenter(test, bary);
  Real first_bary[] = {1./6., 1./6., 1./6.};

  if (!Math::are_vector_equal(3, bary.storage(), first_bary))
    return EXIT_FAILURE;

  test.element = 1;
  interface_mesh.getBarycenter(test, bary);
  Real second_bary[] = {2./3., 2./3., 2./3.};

  if (!Math::are_vector_equal(3, bary.storage(), second_bary))
    return EXIT_FAILURE;

  // Testing a segment completely inside one element
  Point a(0.05, 0.05, 0.05),
        b(0.06, 0.06, 0.06);
  Segment inside_segment(a, b);

  intersector.computeIntersectionQuery(inside_segment);
  test.element = interface_mesh.getNbElement(_segment_2) - 1;
  interface_mesh.getBarycenter(test, bary);
  Real third_bary[] = {0.055, 0.055, 0.055};

  if (!Math::are_vector_equal(3, bary.storage(), third_bary))
    return EXIT_FAILURE;

  // Testing a segment whose end points are inside elements
  Point c(0.1, 0.1, 0.1),
        d(0.9, 0.9, 0.9);
  Segment crossing_segment(c, d);

  intersector.computeIntersectionQuery(crossing_segment);
  UInt el1 = interface_mesh.getNbElement(_segment_2) - 1;
  UInt el2 = el1 - 1;

  test.element = el1;
  interface_mesh.getBarycenter(test, bary);
  Real fourth_bary[] = {13./60., 13./60., 13./60.};

  if (!Math::are_vector_equal(3, bary.storage(), fourth_bary))
    return EXIT_FAILURE;

  test.element = el2;
  interface_mesh.getBarycenter(test, bary);
  Real fifth_bary[] = {37./60., 37./60., 37./60.};

  if (!Math::are_vector_equal(3, bary.storage(), fifth_bary))
    return EXIT_FAILURE;

  // Testing a segment along the edge of elements
  Point e(1, 0, 0),
        f(0, 1, 0);
  Segment edge_segment(e, f);

  UInt current_nb_elements = interface_mesh.getNbElement(_segment_2);

  intersector.computeIntersectionQuery(edge_segment);

  if (interface_mesh.getNbElement(_segment_2) != current_nb_elements + 1)
    return EXIT_FAILURE;

  test.element = interface_mesh.getNbElement(_segment_2) - 1;
  interface_mesh.getBarycenter(test, bary);
  Real sixth_bary[] = {0.5, 0.5, 0};

  if (!Math::are_vector_equal(3, bary.storage(), sixth_bary))
    return EXIT_FAILURE;

  return EXIT_SUCCESS;
}
