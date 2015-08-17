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
add_mesh(swiss_train swiss_train.geo 2 1)

add_library(locomotive_tools
  locomotive_tools.cc
  locomotive_tools.hh
  )

register_example(dumper_low_level
  dumper_low_level.cc
  USE IOHelper
  DEPENDS swiss_train locomotive_tools
  DIRECTORIES_TO_CREATE paraview
  )

register_example(dumpable_interface
  dumpable_interface.cc
  USE IOHelper
  DEPENDS swiss_train locomotive_tools
  DIRECTORIES_TO_CREATE paraview
  )

