/**
 * @file   tetrahedron_tmpl.hh
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Tue Mar 10 2015
 * @date last modification: Tue Mar 10 2015
 *
 * @brief  Tetrahedron class (geometry) for AABB CGAL algos, template impl
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

#ifndef __AKANTU_TETRAHEDRON_TMPL_HH__
#define __AKANTU_TETRAHEDRON_TMPL_HH__

#include "tetrahedron.hh"
#include <CGAL/Cartesian.h>

__BEGIN_AKANTU__

template<typename K>
Tetrahedron<K>::Tetrahedron() :
  CGAL::Tetrahedron_3<K>(),
  meshId(0)
{}

template<typename K>
Tetrahedron<K>::Tetrahedron(const CGAL::Point_3<K> & a,
                            const CGAL::Point_3<K> & b,
                            const CGAL::Point_3<K> & c,
                            const CGAL::Point_3<K> & d) :
  CGAL::Tetrahedron_3<K>(a, b, c, d),
  meshId(0)
{}

template<typename K>
Tetrahedron<K>::Tetrahedron(const Tetrahedron & other) :
  CGAL::Tetrahedron_3<K>(other),
  meshId(other.meshId)
{}

template<typename K>
UInt Tetrahedron<K>::id() const {
  return meshId;
}

template<typename K>
void Tetrahedron<K>::setId(UInt newId) {
  meshId = newId;
}

__END_AKANTU__

#endif
