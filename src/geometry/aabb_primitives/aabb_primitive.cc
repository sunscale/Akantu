/**
 * @file   aabb_primitive.cc
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Clement Roux <clement.roux@epfl.ch>
 *
 * @date creation: Fri Jan 04 2013
 * @date last modification: Wed Jan 31 2018
 *
 * @brief  Macro classe (primitive) for AABB CGAL algos
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

#include "aabb_primitive.hh"

namespace akantu {

Triangle_primitive::Point Triangle_primitive::reference_point() const {
  return primitive.vertex(0);
}

Line_arc_primitive::Point Line_arc_primitive::reference_point() const {
  Real x = to_double(primitive.source().x());
  Real y = to_double(primitive.source().y());
  Real z = to_double(primitive.source().z());
  return cgal::Spherical::Point_3(x, y, z);
}

} // namespace akantu
