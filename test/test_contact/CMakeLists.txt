#===============================================================================
# @file   CMakeLists.txt
#
# @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
#
# @date creation: Tue May 13 2014
# @date last modification: Mon Sep 15 2014
#
# @brief  configuration file for contact tests
#
# @section LICENSE
#
# Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

add_mesh(hertz_2D_mesh_test hertz_2D.geo 2 1)

register_test(test_hertz_2D
  SOURCES hertz_2D.cc
  DEPENDENCIES hertz_2D_mesh_test
  FILES_TO_COPY steel.dat
)

add_mesh(hertz_3D_mesh_test hertz_3D.geo 3 1)

register_test(test_hertz_3D
  SOURCES hertz_3D.cc
  FILES_TO_COPY steel.dat hertz_3D.msh
)

add_mesh(offset_1slave_mesh offset_1slave.geo 2 1)

register_test(test_offset_1slave
  SOURCES offset_1slave.cc
  DEPENDENCIES offset_1slave_mesh
  FILES_TO_COPY steel.dat
)

add_mesh(offset_2slaves_mesh offset_2slaves.geo 2 1)

register_test(test_offset_2slaves
  SOURCES offset_2slaves.cc
  DEPENDENCIES offset_2slaves_mesh
  FILES_TO_COPY steel.dat
)

#===============================================================================
# Alain Curnier suggested tests

add_mesh(acurnier_2D_1_mesh acurnier_2D_1.geo 2 1)

register_test(test_acurnier_2D_1
  SOURCES acurnier_2D_1.cc
  DEPENDENCIES acurnier_2D_1_mesh
  FILES_TO_COPY material.dat
)

add_mesh(acurnier_2D_2_mesh acurnier_2D_2.geo 2 1)

register_test(test_acurnier_2D_2
  SOURCES acurnier_2D_2.cc
  DEPENDENCIES acurnier_2D_2_mesh
  FILES_TO_COPY material.dat
)

add_mesh(acurnier_2D_3_mesh acurnier_2D_3.geo 2 1)

register_test(test_acurnier_2D_3
  SOURCES acurnier_2D_3.cc
  DEPENDENCIES acurnier_2D_3_mesh
  FILES_TO_COPY material.dat
)


add_mesh(acurnier_3D_1_mesh acurnier_3D_1.geo 3 1)

register_test(test_acurnier_3D_1
  SOURCES acurnier_3D_1.cc
  DEPENDENCIES acurnier_3D_1_mesh
  FILES_TO_COPY material.dat
)

add_mesh(acurnier_3D_2_mesh acurnier_3D_2.geo 3 1)

register_test(test_acurnier_3D_2
  SOURCES acurnier_3D_2.cc
  DEPENDENCIES acurnier_3D_2_mesh
  FILES_TO_COPY material.dat
)

add_mesh(acurnier_3D_3_mesh acurnier_3D_3.geo 3 1)

register_test(test_acurnier_3D_3
  SOURCES acurnier_3D_3.cc
  DEPENDENCIES acurnier_3D_3_mesh
  FILES_TO_COPY material.dat
)

