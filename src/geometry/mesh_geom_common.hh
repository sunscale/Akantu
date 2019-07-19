/**
 * @file   mesh_geom_common.hh
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Clement Roux <clement.roux@epfl.ch>
 *
 * @date creation: Fri Jan 04 2013
 * @date last modification: Wed Jan 31 2018
 *
 * @brief  Common file for MeshGeom module
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

#ifndef __AKANTU_MESH_GEOM_COMMON_HH__
#define __AKANTU_MESH_GEOM_COMMON_HH__

#include "aka_common.hh"

#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>

#include <CGAL/Algebraic_kernel_for_spheres_2_3.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Spherical_kernel_3.h>

namespace akantu {

namespace cgal {
  using Cartesian = CGAL::Simple_cartesian<Real>;

  using Spherical = CGAL::Spherical_kernel_3<
      CGAL::Simple_cartesian<CGAL::Quotient<CGAL::MP_Float>>,
      CGAL::Algebraic_kernel_for_spheres_2_3<CGAL::Quotient<CGAL::MP_Float>>>;
} // namespace cgal

} // namespace akantu

#endif // __AKANTU_MESH_GEOM_COMMON_HH__
