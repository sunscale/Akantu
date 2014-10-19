#===============================================================================
# @file   AkantuCPack.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Wed Oct 17 2012
# @date last modification: Tue May 13 2014
#
# @brief  Configure the packaging system
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

set(PACKAGE_FILE_NAME "akantu" CACHE STRING "Name of package to be generated")
mark_as_advanced(PACKAGE_FILE_NAME)

#set(CPACK_GENERATOR "DEB;TGZ;TBZ2;STGZ;RPM")
set(CPACK_GENERATOR "TGZ")

# General configuration
set(CPACK_PACKAGE_VENDOR "LSMS")
set(CPACK_PACKAGE_FILE_NAME "${PACKAGE_FILE_NAME}-${AKANTU_VERSION}-${CPACK_DEBIAN_PACKAGE_ARCHITECTURE}")
set(CPACK_PACKAGE_VERSION "${AKANTU_VERSION}")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "A multipurpose finite element library, Akantu")
set(CPACK_PACKAGE_NAME "akantu")

# Debian config package
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "nicolas.richart@epfl.ch, guillaume.anciaux@epfl.ch")
if(CMAKE_SYSTEM_PROCESSOR MATCHES "i.86" OR CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64")
  if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64" CACHE STRING "Architecture of debian package generation")
  else()
    set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "i386" CACHE STRING "Architecture of debian package generation")
  endif()
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "ppc")
  set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "powerpc" CACHE STRING "Architecture of debian package generation")
else()
  set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "unknown" CACHE STRING "Architecture of debian package generation")
endif()

set(CPACK_DEBIAN_PACKAGE_DEPENDS "${${_project}_PACKAGE_SYSTEM_DEBIAN_PACKAGE_DEPENDS}")
mark_as_advanced(CPACK_DEBIAN_PACKAGE_ARCHITECTURE)
# RPM package configuration
#set(CPACK_RPM_PACKAGE_REQUIRES "${${_project}_PACKAGE_SYSTEM_DEBIAN_PACKAGE_DEPENDS}")


set(CPACK_COMPONENTS_ALL lib dev)
set(CPACK_COMPONENT_LIB_DISPLAY_NAME "Libraries")
set(CPACK_COMPONENT_DEV_DISPLAY_NAME "C++ Headers")
set(CPACK_COMPONENT_DEV_DEPENDS lib)
set(CPACK_COMPONENT_LIB_DESCRIPTION
  "Akantu libraries")
set(CPACK_COMPONENT_DEV_DESCRIPTION
  "Akantu C/C++ header files")
set(CPACK_COMPONENT_LIB_GROUP "Akantu Libraries")
set(CPACK_COMPONENT_DEV_GROUP "Development")

set(CPACK_SOURCE_PACKAGE_FILE_NAME "${PACKAGE_FILE_NAME}-${AKANTU_VERSION}-src")
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/COPYING")

string(TOUPPER ${PROJECT_NAME} _project)

list(APPEND CPACK_SOURCE_IGNORE_FILES ${AKANTU_EXCLUDE_SOURCE_FILES} ${AKANTU_TESTS_EXCLUDE_FILES} ${AKANTU_DOC_EXCLUDE_FILES})
foreach(_pkg ${${_project}_PACKAGE_SYSTEM_PACKAGES_OFF})
  string(TOUPPER "${_pkg}" _pkg)
  list(APPEND CPACK_SOURCE_IGNORE_FILES ${CMAKE_SOURCE_DIR}/packages/${${_project}_${_pkg}_FILE})
endforeach()

list(APPEND CPACK_SOURCE_IGNORE_FILES "/.*build.*/;/CVS/;/\\\\.svn/;/\\\\.bzr/;/\\\\.hg/;/\\\\.hgignore;/\\\\.git/;\\\\.swp$;\\\\.#;/#;~")

include(CPack)
