/**
 * @file   test_geometry_predicates.cc
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Thu Mar 26 2015
 * @date last modification: Thu Mar 26 2015
 *
 * @brief  Tests the geometry predicates
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
#include "geom_helper_functions.hh"

#include <CGAL/Cartesian.h>

#include <iostream>

/* -------------------------------------------------------------------------- */

using namespace akantu;

typedef CGAL::Cartesian<Real> K;
typedef K::Point_3 Point;
typedef K::Segment_3 Segment;

int main(int argc, char * argv[]) {
  debug::setDebugLevel(dblWarning);
  initialize("", argc, argv);

  Point a(0, 1, 0);
  Point b(0, 1, 1);

  Segment seg1(a, b);
  Segment seg2(b, a);

  if (!compareSegments(seg1, seg2))
    return EXIT_FAILURE;

  finalize();
  return EXIT_SUCCESS;
}