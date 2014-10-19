#===============================================================================
# @file   CMakeLists.txt
#
# @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
#
# @date   Mon Mar 24 11:30:00 2014
#
# @brief  configuration file for contact examples
#
# @section LICENSE
#
# Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

register_example(hertz_2D hertz_2D.cc)

add_mesh(hertz_2D_mesh hertz_2D.geo 2 1)
add_dependencies(hertz_2D hertz_2D_mesh)


register_example(hertz_3D hertz_3D.cc)

add_mesh(hertz_3D_mesh hertz_3D.geo 3 1)
add_dependencies(hertz_3D hertz_3D_mesh)


register_example(tetrahedra tetrahedra.cc)

add_mesh(tetrahedra_mesh tetrahedra.geo 3 1)
add_dependencies(tetrahedra tetrahedra_mesh)


register_example(blocks blocks.cc)

add_mesh(blocks_mesh blocks.geo 2 1)
add_dependencies(blocks blocks_mesh)


#===============================================================================
file(COPY steel.dat DESTINATION .)