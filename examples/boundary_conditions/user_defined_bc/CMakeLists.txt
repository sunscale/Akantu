#===============================================================================
# @file   CMakeLists.txt
#
# @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
#
# @date creation: Wed Dec 16 2015
#
# @brief  boundary condition example configuration
#
# @section LICENSE
#
# Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory
# (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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

add_mesh(user_defined_bc_mesh fine_mesh.geo 2 1)

register_example(user_defined_bc
  SOURCES user_defined_bc.cc
  DEPENDS user_defined_bc_mesh
  FILES_TO_COPY material.dat
  )
