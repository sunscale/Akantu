#===============================================================================
# @file   CMakeLists.txt
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Fri Dec 12 2014
# @date last modification: Mon Jan 18 2016
#
# @brief  CMake file for the python wrapping of akantu
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
#===============================================================================

set(PYAKANTU_SRCS
  py_aka_common.cc
  py_aka_error.cc
  py_akantu.cc
  py_boundary_conditions.cc
  py_fe_engine.cc
  py_group_manager.cc
  py_mesh.cc
  py_model.cc
  py_parser.cc
  )

package_is_activated(iohelper _is_activated)
if (_is_activated)
  list(APPEND PYAKANTU_SRCS
    py_dumpable.cc
    )
endif()


package_is_activated(solid_mechanics _is_activated)
if (_is_activated)
  list(APPEND PYAKANTU_SRCS
    py_solid_mechanics_model.cc
    py_material.cc
    )
endif()

package_is_activated(cohesive_element _is_activated)
if (_is_activated)
  list(APPEND PYAKANTU_SRCS
    py_solid_mechanics_model_cohesive.cc
    )
endif()

package_is_activated(heat_transfer _is_activated)
if (_is_activated)
  list(APPEND PYAKANTU_SRCS
    py_heat_transfer_model.cc
    )
endif()

pybind11_add_module(py11_akantu ${PYAKANTU_SRCS})
target_link_libraries(py11_akantu PRIVATE akantu pybind11)
target_include_directories(py11_akantu INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

set_target_properties(py11_akantu PROPERTIES DEBUG_POSTFIX "")
