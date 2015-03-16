/**
 * @file   geom_helper_functions.hh
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Wed Mar 4 2015
 * @date last modification: Thu Mar 5 2015
 *
 * @brief  Helper functions for the computational geometry algorithms
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

#ifndef _AKANTU_GEOM_HELPER_FUNCTIONS_HH__
#define _AKANTU_GEOM_HELPER_FUNCTIONS_HH__

#include "aka_common.hh"

#include <CGAL/Cartesian.h>

__BEGIN_AKANTU__

typedef CGAL::Cartesian<Real> K;

#define EPS 1e-10

inline bool comparePoints(const K::Point_3 & a, const K::Point_3 & b) {
  Math::setTolerance(EPS);
  return Math::are_float_equal(a.x(), b.x()) && Math::are_float_equal(a.y(), b.y());
}

inline bool compareSegments(const K::Segment_3 & a, const K::Segment_3 & b) {
  Math::setTolerance(EPS);
  return (comparePoints(a.source(), b.source()) && comparePoints(a.target(), b.target())) ||
         (comparePoints(a.source(), b.target()) && comparePoints(a.target(), b.source()));
}

struct CompareSegments {
  bool operator()(const K::Segment_3 & a, const K::Segment_3 & b) {
    return compareSegments(a, b);
  }
};

struct CompareSegmentPairs {
  bool operator()(const std::pair<K::Segment_3, UInt> & a, const std::pair<K::Segment_3, UInt> & b) {
    return compareSegments(a.first, b.first);
  }
};
__END_AKANTU__

#endif // _AKANTU_GEOM_HELPER_FUNCTIONS_HH__

