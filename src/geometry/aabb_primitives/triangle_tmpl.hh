/**
 * @file   triangle_tmpl.hh
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Tue Mar 3 2015
 * @date last modification: Tue Mar 3 2015
 *
 * @brief  Triangle classe (geometry) for AABB CGAL algos, template impl
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

#ifndef __AKANTU_TRIANGLE_TMPL_HH__
#define __AKANTU_TRIANGLE_TMPL_HH__

#include "triangle.hh"
#include <CGAL/Cartesian.h>

__BEGIN_AKANTU__

template<typename K>
Triangle<K>::Triangle() :
  CGAL::Triangle_3<K>(),
  meshId(0)
{}

template<typename K>
Triangle<K>::Triangle(const CGAL::Point_3<K> & a, const CGAL::Point_3<K> & b, const CGAL::Point_3<K> & c) : 
  CGAL::Triangle_3<K>(a, b, c),
  meshId(0)
{}

template<typename K>
Triangle<K>::Triangle(const Triangle & other) :
  CGAL::Triangle_3<K>(other),
  meshId(other.meshId)
{}

template<typename K>
UInt Triangle<K>::id() const {
  return meshId;
}

template<typename K>
void Triangle<K>::setId(UInt newId) {
  meshId = newId;
}

__END_AKANTU__

#endif // __AKANTU_TRIANGLE_TMPL_HH__
