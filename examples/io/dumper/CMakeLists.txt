#===============================================================================
# @file   CMakeLists.txt
#
# @author Fabian Barras <fabian.barras@epfl.ch>
#
# @date creation: Fri Sep 03 2010
# @date last modification: Wed Jan 20 2016
#
# @brief  CMakeLists for DumperIOHelper examples
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

add_mesh(swiss_train_mesh swiss_train.geo 2 1)

add_library(locomotive_tools
  locomotive_tools.cc
  locomotive_tools.hh
  )

package_get_include_dir(BOOST _boost_include_dir)

target_link_libraries(locomotive_tools PRIVATE akantu)
target_include_directories(locomotive_tools PRIVATE ${AKANTU_INCLUDE_DIRS} ${_boost_include_dir})

if(AKANTU_EXTRA_CXX_FLAGS)
  set_target_properties(locomotive_tools PROPERTIES COMPILE_FLAGS ${AKANTU_EXTRA_CXX_FLAGS})
endif()

register_example(dumper_low_level
  SOURCES dumper_low_level.cc
  USE_PACKAGES IOHelper
  DEPENDS swiss_train_mesh locomotive_tools
  DIRECTORIES_TO_CREATE paraview
  )

register_example(dumpable_interface
  SOURCES dumpable_interface.cc
  USE_PACKAGES IOHelper
  DEPENDS swiss_train_mesh locomotive_tools
  DIRECTORIES_TO_CREATE paraview
  )

