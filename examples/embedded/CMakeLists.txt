#===============================================================================
# @file   CMakeLists.txt
#
# @author Lucas Frérot <lucas.frerot@epfl.ch>
#
# @date creation: Wed Jul 22 2015
# @date last modification: Wed Jul 22 2015
#
# @brief  configuration for embedded example
#
# @section LICENSE
#
# Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# Akantu is free  software: you can redistribute it and/or  modify it under the
# terms  of the  GNU Lesser  General Public  License as  published by  the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
# details.
#
# You should  have received  a copy  of the GNU  Lesser General  Public License
# along with Akantu. If not, see <http://www.gnu.org/licenses/>.
#
# @section DESCRIPTION
#
#===============================================================================

register_example(embedded embedded.cc USE CGAL)

add_mesh(concrete_mesh concrete.geo 2 1)
add_mesh(reinforcement_mesh reinforcement.geo 1 1)
add_dependencies(embedded concrete_mesh)
add_dependencies(embedded reinforcement_mesh)

#===============================================================================
file(COPY material.dat DESTINATION .)

