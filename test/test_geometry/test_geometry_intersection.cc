/**
 * @file   test_geometry_intersection.cc
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Fri Feb 27 2015
 * @date last modification: Thu Mar 5 2015
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
#include "mesh_geom_factory.hh"
#include "tree_type_helper.hh"
#include "geom_helper_functions.hh"

#include "mesh_geom_common.hh"

#include <iostream>

/* -------------------------------------------------------------------------- */

using namespace akantu;

typedef Cartesian K;
typedef IntersectionTypeHelper<TreeTypeHelper<Triangle<K>, K>, K::Segment_3>::intersection_type result_type;

/* -------------------------------------------------------------------------- */

int main (int argc, char * argv[]) {
  initialize("", argc, argv);
  debug::setDebugLevel(dblWarning);

  Mesh mesh(2);
  mesh.read("test_geometry_triangle.msh");

  MeshGeomFactory<2, _triangle_3, Triangle<K>, K> factory(mesh);
  factory.constructData();

  const TreeTypeHelper<Triangle<K>, K>::tree & tree = factory.getTree();

  K::Point_3 a(0., 0.25, 0.), b(1., 0.25, 0.);
  K::Segment_3 line(a, b);

  K::Point_3 begin(a), intermediate(0.25, 0.25, 0.), end(0.75, 0.25, 0.);
  K::Segment_3 result_0(begin, intermediate), result_1(intermediate, end);

  std::list<result_type> list_of_intersections;
  tree.all_intersections(line, std::back_inserter(list_of_intersections));

  const result_type & intersection_0 = list_of_intersections.back();
  const result_type & intersection_1 = list_of_intersections.front();

  if (!intersection_0 || !intersection_1)
    return EXIT_FAILURE;

  /// *-> first is the intersection ; *->second is the primitive id
  if (const K::Segment_3 * segment = boost::get<K::Segment_3>(&(intersection_0->first))) {
    if (!compareSegments(*segment, result_0)) {
      return EXIT_FAILURE;
    }
  } else return EXIT_FAILURE;

  if (const K::Segment_3 * segment = boost::get<K::Segment_3>(&(intersection_1->first))) {
    if (!compareSegments(*segment, result_1)) {
      return EXIT_FAILURE;
    }
  } else return EXIT_FAILURE;

  finalize();
  return EXIT_SUCCESS;
}

