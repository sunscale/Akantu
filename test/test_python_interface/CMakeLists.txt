#===============================================================================
# @file   CMakeLists.txt
#
# @author Fabian Barras <fabian.barras@epfl.ch>
# @author Lucas Frerot <lucas.frerot@epfl.ch>
#
# @date creation: Fri Sep 03 2010
# @date last modification: Mon Feb 05 2018
#
# @brief  Python Interface tests
#
# @section LICENSE
#
# Copyright (©) 2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# Akantu is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with Akantu. If not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================
#get_target_property(_srcs py11_akantu SOURCES)
#foreach(_src ${_srcs})
#  list(APPEND py11_akantu_test_common_SRCS ${PROJECT_SOURCE_DIR}/python/${_src})
#endforeach()
#
#akantu_pybind11_add_module(py11_akantu_test_common MODULE
#  test_common.cc
#  ${py11_akantu_test_common_SRCS})
#target_link_libraries(py11_akantu_test_common
#  PRIVATE $<TARGET_PROPERTY:py11_akantu,LINK_LIBRARIES>)
#target_include_directories(py11_akantu_test_common
#  PRIVATE $<TARGET_PROPERTY:py11_akantu,INCLUDE_DIRECTORIES>)
#

akantu_pybind11_add_module(py11_akantu_test_common MODULE test_common.cc)

add_mesh(mesh_dcb_2d mesh_dcb_2d.geo 2 2)

register_test(test_python_interface
  SCRIPT test_pybind.py
  PYTHON
  FILES_TO_COPY elastic.dat
  DEPENDS mesh_dcb_2d py11_akantu_test_common
  PACKAGE python_interface
  )
