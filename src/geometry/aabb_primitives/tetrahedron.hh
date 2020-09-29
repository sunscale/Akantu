/**
 * @file   tetrahedron.hh
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 *
 * @date creation: Fri Feb 27 2015
 * @date last modification: Mon Jun 19 2017
 *
 * @brief  Tetrahedron classe (geometry) for AABB CGAL algos
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

#ifndef AKANTU_TETRAHEDRON_HH_
#define AKANTU_TETRAHEDRON_HH_

#include "aka_common.hh"

#include "mesh_geom_common.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */

/// Class used for substitution of CGAL::Tetrahedron_3 primitive
template <typename K> class Tetrahedron : public CGAL::Tetrahedron_3<K> {
public:
  /// Default constructor
  Tetrahedron() : CGAL::Tetrahedron_3<K>() {}

  /// Copy constructor
  Tetrahedron(const Tetrahedron & other)
      : CGAL::Tetrahedron_3<K>(other), meshId(other.meshId) {}

  /// Construct from 4 points
  Tetrahedron(const CGAL::Point_3<K> & a, const CGAL::Point_3<K> & b,
              const CGAL::Point_3<K> & c, const CGAL::Point_3<K> & d)
      : CGAL::Tetrahedron_3<K>(a, b, c, d) {}

public:
  UInt id() const { return meshId; }
  void setId(UInt newId) { meshId = newId; }

protected:
  /// Id of the element represented by the primitive
  UInt meshId{0};
};

} // namespace akantu

#endif
