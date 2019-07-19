/**
 * @file   geom_helper_functions.hh
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Clement Roux <clement.roux@epfl.ch>
 *
 * @date creation: Fri Jan 04 2013
 * @date last modification: Wed Jan 31 2018
 *
 * @brief  Helper functions for the computational geometry algorithms
 *
 * @section LICENSE
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

#ifndef _AKANTU_GEOM_HELPER_FUNCTIONS_HH__
#define _AKANTU_GEOM_HELPER_FUNCTIONS_HH__

#include "aka_common.hh"
#include "aka_math.hh"
#include "tree_type_helper.hh"

#include "mesh_geom_common.hh"

namespace akantu {

/// Fuzzy compare of two points
template <class Point>
inline bool comparePoints(const Point & a, const Point & b) {
  return Math::are_float_equal(a.x(), b.x()) &&
         Math::are_float_equal(a.y(), b.y()) &&
         Math::are_float_equal(a.z(), b.z());
}

template <>
inline bool comparePoints(const cgal::Spherical::Circular_arc_point_3 & a,
                          const cgal::Spherical::Circular_arc_point_3 & b) {
  return Math::are_float_equal(CGAL::to_double(a.x()),
                               CGAL::to_double(b.x())) &&
         Math::are_float_equal(CGAL::to_double(a.y()),
                               CGAL::to_double(b.y())) &&
         Math::are_float_equal(CGAL::to_double(a.z()), CGAL::to_double(b.z()));
}

/// Fuzzy compare of two segments
template <class K>
inline bool compareSegments(const CGAL::Segment_3<K> & a,
                            const CGAL::Segment_3<K> & b) {
  return (comparePoints(a.source(), b.source()) &&
          comparePoints(a.target(), b.target())) ||
         (comparePoints(a.source(), b.target()) &&
          comparePoints(a.target(), b.source()));
}

/// Compare segment pairs
inline bool
compareSegmentPairs(const std::pair<cgal::Cartesian::Segment_3, UInt> & a,
                    const std::pair<cgal::Cartesian::Segment_3, UInt> & b) {
  return compareSegments(a.first, b.first);
}

/// Pair ordering operator based on first member
struct segmentPairsLess {
  inline bool
  operator()(const std::pair<cgal::Cartesian::Segment_3, UInt> & a,
             const std::pair<cgal::Cartesian::Segment_3, UInt> & b) {
    return CGAL::compare_lexicographically(a.first.min(), b.first.min()) ||
           CGAL::compare_lexicographically(a.first.max(), b.first.max());
  }
};

/* -------------------------------------------------------------------------- */
/* Predicates                                                                 */
/* -------------------------------------------------------------------------- */

/// Predicate used to determine if two segments are equal
class IsSameSegment {

public:
  IsSameSegment(const cgal::Cartesian::Segment_3 & segment)
      : segment(segment) {}

  bool
  operator()(const std::pair<cgal::Cartesian::Segment_3, UInt> & test_pair) {
    return compareSegments(segment, test_pair.first);
  }

protected:
  const cgal::Cartesian::Segment_3 segment;
};

} // namespace akantu

#endif // _AKANTU_GEOM_HELPER_FUNCTIONS_HH__
