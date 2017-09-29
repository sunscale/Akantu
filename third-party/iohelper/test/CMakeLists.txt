#===============================================================================
# @file   CMakeLists.txt
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
#
# @date creation: Thu Oct 11 2012
# @date last modification: Tue May 14 2013
#
# @brief  configuration file for tests
#
# @section LICENSE
#
# Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# IOHelper is free  software: you can redistribute it and/or  modify it under the
# terms  of the  GNU Lesser  General Public  License as  published by  the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# IOHelper is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
# details.
#
# You should  have received  a copy  of the GNU  Lesser General  Public License
# along with IOHelper. If not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================

SET(test_paraview_dumper_SRC
  mesh_io_msh.cc
  test_paraview_dumper.cc
  )

SET(test_text_dumper_SRC
  mesh_io_msh.cc
  test_text_dumper.cc
  )

add_executable(test_paraview_dumper ${test_paraview_dumper_SRC})
add_executable(test_text_dumper ${test_text_dumper_SRC})
target_link_libraries(test_paraview_dumper iohelper)
target_link_libraries(test_text_dumper iohelper)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/cube.msh DESTINATION .)

add_test(test_paraview_dumper test_paraview_dumper cube.msh 3)
add_test(test_text_dumper test_text_dumper cube.msh 3)

