/**
 * @file   test_geometry_intersection.cc
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Fri Feb 27 2015
 * @date last modification: Fri Feb 27 2015
 *
 * @brief  Tests the intersection module
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
#include "mesh_tree_constructor.hh"
#include "tree_type_helper.hh"

#include <CGAL/Cartesian.h>

#include <iostream>

using namespace akantu;

typedef CGAL::Cartesian<Real> K;

int main (int argc, char * argv[]) {
  debug::setDebugLevel(dblWarning);
  initialize("", argc, argv);

  Mesh mesh(2);
  mesh.read("mesh.msh");

  MeshTreeConstructor<2, _triangle_3> tree_constructor(mesh);
  tree_constructor.constructData();

  const TreeTypeHelper<2, _triangle_3>::tree & tree = tree_constructor.getTree();

  K::Point_3 a(0., 0.25, 0.), b(1., 0.25, 0.);
  K::Line_3 line(a, b);

  std::cout << line << std::endl;
  std::cout << tree.number_of_intersected_primitives(line) << std::endl;

  boost::optional < TreeTypeHelper<2, _triangle_3>::tree::Intersection_and_primitive_id<K::Line_3>::Type > segment_intersection = tree.any_intersection(line);

  if (segment_intersection) {
    if (const K::Segment_3 * segment = boost::get<K::Segment_3>(&(segment_intersection->first))) {
      std::cout << *segment << std::endl;
    }
  }

  finalize();
  return 0;
}
