#===============================================================================
# @file   CMakeLists.txt
#
# @author Fabian Barras <fabian.barras@epfl.ch>
#
# @date creation: Fri Aug 14 2015
#
# @brief  CMakeLists for DumperIOHelper examples
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

register_example(example_dumper_low_level
  SOURCES example_dumper_low_level.cc
  FILES_TO_COPY swiss_train.msh
  DIRECTORIES_TO_CREATE paraview
  PACKAGE iohelper
  )

register_example(example_dumpable_interface
  SOURCES example_dumpable_interface.cc
  FILES_TO_COPY swiss_train.msh
  DIRECTORIES_TO_CREATE paraview
  PACKAGE iohelper
  )

#===============================================================================
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview)