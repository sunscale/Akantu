/**
 * @file mesh_geom_common.hh
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 *
 * @date creation: Wed May 13 2015
 * @date last modification: Wed May 13 2015
 *
 * @brief Common file for MeshGeom module
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2015 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef __AKANTU_MESH_GEOM_COMMON_HH__
#define __AKANTU_MESH_GEOM_COMMON_HH__

#include "aka_common.hh"

#include <CGAL/Cartesian.h>
#include <CGAL/Simple_cartesian.h>

__BEGIN_AKANTU__

typedef CGAL::Simple_cartesian<double> Cartesian;

__END_AKANTU__

#endif // __AKANTU_MESH_GEOM_COMMON_HH__
