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
  debug::setDebugLevel(dblWarning);
  initialize("", argc, argv);

  Mesh mesh(3), interface_mesh(3, "interface_mesh");
  mesh.read("test_geometry_tetrahedron.msh");

  MeshSegmentIntersector<3, _tetrahedron_4> intersector(mesh, interface_mesh);
  intersector.constructData();

  Point point(1., 1., 1.);
  Segment segment(CGAL::ORIGIN, point);

  intersector.computeIntersectionQuery(segment);

  if (interface_mesh.getNbElement(_segment_2) != 4)
    return EXIT_FAILURE;

  return EXIT_SUCCESS;
}
