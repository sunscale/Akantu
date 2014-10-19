/**
 * @file   aka_geometry.cc
 *
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 *
 * @date creation: Fri Jan 04 2013
 * @date last modification: Tue May 13 2014
 *
 * @brief  geometric operations
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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


#include "aka_geometry.hh"

__BEGIN_AKANTU__

using std::cout;
using std::endl;

/*! \param p - A constant reference to the first point that defines the line segment.
 * \param q - A constant reference to the second point that defines the line segment.
 * \param r - A constant reference to the point for which the predicate is computed.
 * \return A double, being -1.0 if the point lines on one side of the line and 1.0 if it lies
 * on the other side.
 */
Real left_turn(const Point<2>& p, const Point<2>& q, const Point<2>& r) {
    
  if(((q[0]-p[0]) * (r[1]-p[1])) > ((r[0]-p[0]) * (q[1]-p[1])))
    return 1.;
  else
    return -1.;
}


__END_AKANTU__
