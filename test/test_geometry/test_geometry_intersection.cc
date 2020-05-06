/**
 * @file   test_geometry_intersection.cc
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Clement Roux <clement.roux@epfl.ch>
 *
 * @date creation: Fri Feb 27 2015
 * @date last modification: Wed Jan 31 2018
 *
 * @brief  Tests the intersection module
 *
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "geom_helper_functions.hh"
#include "mesh_geom_factory.hh"
#include "tree_type_helper.hh"

#include "mesh_geom_common.hh"

#include <iostream>
#include <iterator>

#include <CGAL/Exact_spherical_kernel_3.h>
#include <CGAL/Spherical_kernel_intersections.h>
#include <CGAL/intersections.h>

/* -------------------------------------------------------------------------- */

using namespace akantu;

typedef cgal::Cartesian K;
typedef IntersectionTypeHelper<TreeTypeHelper<Triangle<K>, K>,
                               K::Segment_3>::intersection_type result_type;

typedef cgal::Spherical SK;
typedef boost::variant<std::pair<SK::Circular_arc_point_3, UInt>> sk_inter_res;
/*typedef CGAL::cpp11::result_of<SK::Intersect_3(SK::Line_arc_3,
                           SK::Sphere_3,
                           std::back_insert_iterator<
                           std::list<sk_inter_res> >)>::type sk_res;*/

typedef std::pair<SK::Circular_arc_point_3, UInt> pair_type;

/* -------------------------------------------------------------------------- */

int main(int argc, char * argv[]) {
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
  if (const K::Segment_3 * segment =
          boost::get<K::Segment_3>(&(intersection_0->first))) {
    if (!compareSegments(*segment, result_0)) {
      return EXIT_FAILURE;
    }
  } else
    return EXIT_FAILURE;

  if (const K::Segment_3 * segment =
          boost::get<K::Segment_3>(&(intersection_1->first))) {
    if (!compareSegments(*segment, result_1)) {
      return EXIT_FAILURE;
    }
  } else
    return EXIT_FAILURE;

  SK::Sphere_3 sphere(SK::Point_3(0, 0, 0), 3.);
  SK::Segment_3 seg(SK::Point_3(0, 0, 0), SK::Point_3(2., 2., 2.));
  SK::Line_arc_3 arc(seg);

  std::list<sk_inter_res> s_results;
  CGAL::intersection(arc, sphere, std::back_inserter(s_results));

  if (pair_type * pair = boost::get<pair_type>(&s_results.front())) {
    std::cout << "xi = " << to_double(pair->first.x())
              << ", yi = " << to_double(pair->first.y()) << std::endl;
    if (!comparePoints(pair->first, SK::Circular_arc_point_3(1.0, 1.0, 1.0)))
      return EXIT_FAILURE;
  } else
    return EXIT_FAILURE;

  MeshGeomFactory<2, _triangle_3, Line_arc<SK>, SK> Sfactory(mesh);
  Sfactory.constructData();

  finalize();
  return EXIT_SUCCESS;
}
