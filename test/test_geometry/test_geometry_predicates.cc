/**
 * @file   test_geometry_predicates.cc
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 *
 * @date creation: Fri Jan 04 2013
 * @date last modification: Wed Jan 31 2018
 *
 * @brief  Tests the geometry predicates
 *
 *
 * Copyright (©) 2014-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#include "mesh_geom_common.hh"

#include <iostream>

/* -------------------------------------------------------------------------- */

using namespace akantu;

typedef cgal::Cartesian K;
typedef K::Point_3 Point;
typedef K::Segment_3 Segment;

int main(int argc, char * argv[]) {
  initialize("", argc, argv);
  debug::setDebugLevel(dblWarning);

  Point a(0, 1, 0);
  Point b(0, 1, 1);

  Segment seg1(a, b);
  Segment seg2(b, a);

  if (!compareSegments(seg1, seg2))
    return EXIT_FAILURE;

  // Testing sort + unique on list of segments
  std::vector<std::pair<K::Segment_3, UInt>> pair_list;
  pair_list.push_back(std::make_pair(seg1, 1));
  pair_list.push_back(std::make_pair(seg2, 2));

  segmentPairsLess sorter;
  std::sort(pair_list.begin(), pair_list.end(), sorter);
  std::vector<std::pair<K::Segment_3, UInt>>::iterator it =
      std::unique(pair_list.begin(), pair_list.end(), compareSegmentPairs);

  if (it - pair_list.begin() != 1) {
    std::cout << pair_list.size() << std::endl;
    return EXIT_FAILURE;
  }

  // Testing insertion in set
  std::set<std::pair<K::Segment_3, UInt>, segmentPairsLess> pair_set;
  pair_set.insert(pair_set.begin(), std::make_pair(seg1, 1));
  pair_set.insert(pair_set.begin(), std::make_pair(seg2, 2));

  if (pair_set.size() != 1) {
    std::cout << pair_set.size() << std::endl;
    return EXIT_FAILURE;
  }

  finalize();
  return EXIT_SUCCESS;
}
