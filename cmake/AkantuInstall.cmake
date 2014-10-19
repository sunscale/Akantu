#===============================================================================
# @file   AkantuInstall.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Wed Oct 17 2012
# @date last modification: Fri Sep 19 2014
#
# @brief  Create the files that allows users to link with Akantu in an other
# cmake project
#
# @section LICENSE
#
# Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#===============================================================================

#===============================================================================
# Config gen for external packages
#===============================================================================
configure_file(cmake/AkantuBuildTreeSettings.cmake.in  "${PROJECT_BINARY_DIR}/AkantuBuildTreeSettings.cmake" @ONLY)

file(WRITE "${PROJECT_BINARY_DIR}/AkantuConfigInclude.cmake" "
#===============================================================================
# @file   AkantuConfigInclude.cmake
# @author Nicolas Richart <nicolas.richart@epfl.ch>
# @date   Fri Jun 11 09:46:59 2010
#
# @section LICENSE
#
# Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

")

foreach(_option ${AKANTU_PACKAGE_SYSTEM_PACKAGES_NAMES_LIST_ALL})
  list(FIND AKANTU_OPTION_LIST ${_option_name} _index)
  if (_index EQUAL -1)
    if(NOT "${_option}" STREQUAL "CORE")
      if(NOT AKANTU_${_option})
        set(AKANTU_${_option} OFF)
      endif()
      file(APPEND "${PROJECT_BINARY_DIR}/AkantuConfigInclude.cmake" "
set(AKANTU_HAS_${_option} ${AKANTU_${_option}})")
    endif()
  endif()
endforeach()

file(APPEND "${PROJECT_BINARY_DIR}/AkantuConfigInclude.cmake"
"

set(AKANTU_HAS_PARTITIONER  ${AKANTU_PARTITIONER})
set(AKANTU_HAS_SOLVER       ${AKANTU_SOLVER})
")

foreach(_option ${AKANTU_OPTION_LIST})
  package_pkg_name(${_option} _pkg_name)
  file(APPEND "${PROJECT_BINARY_DIR}/AkantuConfigInclude.cmake" "
list(APPEND AKANTU_OPTION_LIST ${_option})
set(AKANTU_USE_${_option} ${AKANTU_${_option}})")
  if(${_pkg_name}_LIBRARIES)
    file(APPEND "${PROJECT_BINARY_DIR}/AkantuConfigInclude.cmake" "
set(${_pkg_name}_LIBRARIES ${${_pkg_name}_LIBRARIES})")
  endif()
  if(${_pkg_name}_INCLUDE_DIR)
    file(APPEND "${PROJECT_BINARY_DIR}/AkantuConfigInclude.cmake" "
set(${_pkg_name}_INCLUDE_DIR ${${_pkg_name}_INCLUDE_DIR})
")
  endif()
endforeach()

file(APPEND "${PROJECT_BINARY_DIR}/AkantuConfigInclude.cmake" "
set(AKANTU_BOOST_INCLUDE_DIR ${Boost_INCLUDE_DIRS})
set(AKANTU_BOOST_LIBRARIES ${Boost_LIBRARIES})
")


# Create the AkantuConfig.cmake and AkantuConfigVersion files
get_filename_component(CONF_REL_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}" ABSOLUTE)
configure_file(cmake/AkantuConfig.cmake.in "${PROJECT_BINARY_DIR}/AkantuConfig.cmake" @ONLY)
configure_file(cmake/AkantuConfigVersion.cmake.in "${PROJECT_BINARY_DIR}/AkantuConfigVersion.cmake" @ONLY)
configure_file(cmake/AkantuUse.cmake "${PROJECT_BINARY_DIR}/AkantuUse.cmake" COPYONLY)


# include(CMakePackageConfigHelpers)

# configure_package_config_file(cmake/AkantuConfig.cmake.in ${PROJECT_BINARY_DIR}/AkantuConfig.cmake
#   INSTALL_DESTINATION "${CONF_REL_INCLUDE_DIR}lib/akantu/cmake"
#   PATH_VARS "${CONF_REL_INCLUDE_DIR}/include" )

# write_basic_package_version_file(${PROJECT_BINARY_DIR}/AkantuConfigVersion.cmake
#   VERSION "${AKANTU_VERSION}"
#   COMPATIBILITY SameMajorVersion)


# Install the export set for use with the install-tree
install(FILES ${PROJECT_BINARY_DIR}/AkantuConfig.cmake
  ${PROJECT_BINARY_DIR}/AkantuConfigInclude.cmake
  ${PROJECT_BINARY_DIR}/AkantuConfigVersion.cmake
  ${PROJECT_SOURCE_DIR}/cmake/AkantuUse.cmake
  DESTINATION  lib/akantu
  COMPONENT dev)

install(FILES
  ${PROJECT_SOURCE_DIR}/cmake/Modules/FindIOHelper.cmake
  ${PROJECT_SOURCE_DIR}/cmake/Modules/FindQVIEW.cmake
  ${PROJECT_SOURCE_DIR}/cmake/Modules/FindMumps.cmake
  ${PROJECT_SOURCE_DIR}/cmake/Modules/FindScotch.cmake
  ${PROJECT_SOURCE_DIR}/cmake/Modules/FindGMSH.cmake
  DESTINATION  lib/akantu/cmake
  COMPONENT dev)

