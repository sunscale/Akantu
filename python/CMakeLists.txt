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


add_library(pyakantu py_aka_model.cc py_aka_error.cc py_aka_solid_mechanics_model.cc py_aka_heat_transfer_model.cc py_aka_solid_mechanics_model_cohesive.cc py_aka_material.cc py_aka_boundary_conditions.cc py_aka_common.cc py_aka_mesh.cc py_aka_parser.cc py_aka_fe_engine.cc)
target_link_libraries(pyakantu PUBLIC  akantu pybind11)
pybind11_add_module(py11_akantu py_akantu.cc)
target_link_libraries(py11_akantu PRIVATE pyakantu)

set_property(TARGET py11_akantu PROPERTY DEBUG_POSTFIX "")
