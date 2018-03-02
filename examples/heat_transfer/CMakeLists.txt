#===============================================================================
# @file   CMakeLists.txt
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
#
# @date creation: Mon Jan 18 2016
#
# @brief  configuration for heat transfer example
#
# @section LICENSE
#
# Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
# Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
# Solides)
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

add_mesh(heat_transfer_mesh_2d square.geo 2 2)
add_mesh(heat_transfer_mesh_3d cube.geo 3 2)

#add_mesh(heat_transfer_mesh_2d_struct square.geo 2 2)
#add_mesh(heat_transfer_mesh_3d_struct cube.geo 3 3)
#
register_example(heat_transfer_dynamics_2d
  SOURCES heat_transfer_dynamics_2d.cc
  DEPENDS heat_transfer_mesh_2d
  FILES_TO_COPY material.dat
  )

register_example(heat_transfer_static_2d
  SOURCES heat_transfer_static_2d.cc
  DEPENDS heat_transfer_mesh_2d
  FILES_TO_COPY material.dat
  )

register_example(heat_transfer_dynamics_3d
  SOURCES heat_transfer_dynamics_3d.cc
  DEPENDS heat_transfer_mesh_3d
  FILES_TO_COPY material.dat
  )

