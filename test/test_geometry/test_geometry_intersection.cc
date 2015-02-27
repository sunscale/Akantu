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

/* -------------------------------------------------------------------------- */

using namespace akantu;

typedef CGAL::Cartesian<Real> K;
typedef boost::optional < TreeTypeHelper<2, _triangle_3>::tree::Intersection_and_primitive_id<K::Line_3>::Type > Line_intersection;

/* -------------------------------------------------------------------------- */

bool comparePoints(const K::Point_3 & a, const K::Point_3 & b);
bool compareSegments(const K::Segment_3 & a, const K::Segment_3 & b);

/* -------------------------------------------------------------------------- */

int main (int argc, char * argv[]) {
  debug::setDebugLevel(dblWarning);
  initialize("", argc, argv);

  Math::setTolerance(1e-10);

  Mesh mesh(2);
  mesh.read("mesh.msh");

  MeshTreeConstructor<2, _triangle_3> tree_constructor(mesh);
  tree_constructor.constructData();

  const TreeTypeHelper<2, _triangle_3>::tree & tree = tree_constructor.getTree();

  K::Point_3 a(0., 0.25, 0.), b(1., 0.25, 0.);
  K::Line_3 line(a, b);

  if (tree.number_of_intersected_primitives(line) != 2)
    return EXIT_FAILURE;

  K::Point_3 begin(a), intermediate(0.25, 0.25, 0.), end(0.75, 0.25, 0.);
  K::Segment_3 result_0(begin, intermediate), result_1(intermediate, end);

  std::list<Line_intersection> list_of_intersections;
  tree.all_intersections(line, std::back_inserter(list_of_intersections));

  const Line_intersection & intersection_0 = list_of_intersections.front();
  const Line_intersection & intersection_1 = list_of_intersections.back();

  if (!intersection_0 || !intersection_1)
    return EXIT_FAILURE;

  if (const K::Segment_3 * segment = boost::get<K::Segment_3>(&(intersection_0->first))) {
    if (!compareSegments(*segment, result_0)) {
      return EXIT_FAILURE;
    }
  }

  if (const K::Segment_3 * segment = boost::get<K::Segment_3>(&(intersection_1->first))) {
    if (!compareSegments(*segment, result_1)) {
      return EXIT_FAILURE;
    }
  }

  finalize();
  return EXIT_SUCCESS;
}

/* -------------------------------------------------------------------------- */


bool comparePoints(const K::Point_3 & a, const K::Point_3 & b) {
  return Math::are_float_equal(a.x(), b.x()) && Math::are_float_equal(a.y(), b.y());
}

bool compareSegments(const K::Segment_3 & a, const K::Segment_3 & b) {
  return (comparePoints(a.source(), b.source()) && comparePoints(a.target(), b.target())) ||
         (comparePoints(a.source(), b.target()) && comparePoints(a.target(), b.source()));
}

