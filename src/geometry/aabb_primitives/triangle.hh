/**
 * @file   triangle.hh
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 *
 * @date creation: Fri Jan 04 2013
 * @date last modification: Mon Jun 19 2017
 *
 * @brief  Triangle classe (geometry) for AABB CGAL algos
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

#ifndef __AKANTU_TRIANGLE_HH__
#define __AKANTU_TRIANGLE_HH__

#include "aka_common.hh"

#include "mesh_geom_common.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */

/// Class used for substitution of CGAL::Triangle_3 primitive
template <typename K> class Triangle : public CGAL::Triangle_3<K> {
public:
  /// Default constructor
  Triangle() : CGAL::Triangle_3<K>(), meshId(0) {}

  /// Copy constructor
  Triangle(const Triangle & other)
      : CGAL::Triangle_3<K>(other), meshId(other.meshId) {}

  /// Construct from 3 points
  Triangle(const CGAL::Point_3<K> & a, const CGAL::Point_3<K> & b,
           const CGAL::Point_3<K> & c)
      : CGAL::Triangle_3<K>(a, b, c), meshId(0) {}

public:
  UInt id() const { return meshId; }
  void setId(UInt newId) { meshId = newId; }

protected:
  /// Id of the element represented by the primitive
  UInt meshId;
};

} // namespace akantu

#endif // __AKANTU_TRIANGLE_HH__
