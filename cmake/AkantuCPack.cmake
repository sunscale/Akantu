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
if(NOT CMAKE_SYSTEM_NAME STREQUAL "Windows")
  set(CPACK_GENERATOR "TGZ")
else()
  set(CPACK_GENERATOR "TGZ;NSIS")

  package_get_all_external_informations(
    _external_include_dirs
    _external_libraries
  )
  set(CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS ${_external_libraries})

  include(InstallRequiredSystemLibraries)
endif()

if(CMAKE_SYSTEM_PROCESSOR MATCHES "i.86" OR CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64" OR CMAKE_SYSTEM_PROCESSOR MATCHES "[aA][mM][dD]64")
  if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(_arch "amd64")
  else()
    set(_arch "i386")
  endif()
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "ppc")
  set(_arch "powerpc")
else()
  set(_arch "unknown")
endif()

if(WIN32 AND MINGW)
  set(_arch "mingw32")
endif()

# General configuration
set(CPACK_PACKAGE_VENDOR "LSMS")
set(CPACK_PACKAGE_FILE_NAME "${PACKAGE_FILE_NAME}-${AKANTU_VERSION}-${_arch}")
set(CPACK_PACKAGE_VERSION "${AKANTU_VERSION}")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "A multipurpose finite element library, Akantu")
set(CPACK_PACKAGE_NAME "akantu")
#set(CMAKE_PACKAGE_ICON "${PROJECT_SOURCE_DIR}/cmake/akantu.ico")

# Debian config package
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "nicolas.richart@epfl.ch, guillaume.anciaux@epfl.ch")
set(CPACK_DEBIAN_PACKAGE_SECTION "Science")
set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "${_arch}" CACHE STRING "Architecture of akantu's package")
set(CPACK_DEBIAN_PACKAGE_DEPENDS "${${_project}_PACKAGE_SYSTEM_DEBIAN_PACKAGE_DEPENDS}")
mark_as_advanced(CPACK_DEBIAN_PACKAGE_ARCHITECTURE)
# RPM package configuration
#set(CPACK_RPM_PACKAGE_REQUIRES "${${_project}_PACKAGE_SYSTEM_DEBIAN_PACKAGE_DEPENDS}")

# NSIS Windows installer
#set(CPACK_NSIS_MUI_ICON "${PROJECT_SOURCE_DIR}/cmake/akantu.ico")
#set(CPACK_NSIS_CONTACT "akantu@akantu.ch")
#set(CPACK_NSIS_MODIFY_PATH ON)

# Components description
set(CPACK_COMPONENTS_ALL lib dev bin python)
set(CPACK_COMPONENT_LIB_DISPLAY_NAME "Libraries")
set(CPACK_COMPONENT_BIN_DISPLAY_NAME "Examples")
set(CPACK_COMPONENT_PYTHON_DISPLAY_NAME "Python interface")
set(CPACK_COMPONENT_DEV_DISPLAY_NAME "C++ Headers")
set(CPACK_COMPONENT_DEV_DEPENDS lib)
set(CPACK_COMPONENT_BIN_DEPENDS lib)
set(CPACK_COMPONENT_PYTHON_DEPENDS lib)
set(CPACK_COMPONENT_LIB_DESCRIPTION
  "Akantu libraries")
set(CPACK_COMPONENT_DEV_DESCRIPTION
  "Akantu C/C++ header files")
set(CPACK_COMPONENT_LIB_GROUP "Akantu Libraries and Executables")
set(CPACK_COMPONENT_BIN_GROUP "Akantu Libraries and Executables")
set(CPACK_COMPONENT_PYTHON_GROUP "Akantu Libraries and Executables")
set(CPACK_COMPONENT_DEV_GROUP "Development")

set(CPACK_SOURCE_PACKAGE_FILE_NAME "${PACKAGE_FILE_NAME}-${AKANTU_VERSION}-src")
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/COPYING")

string(TOUPPER ${PROJECT_NAME} _project)

unset(CPACK_SOURCE_IGNORE_FILES)
package_get_all_deactivated_packages(_deactivated_packages)
foreach(_pkg ${_deactivated_packages}})
  _package_get_filename(${_pkg} _file_name)
  list(APPEND CPACK_SOURCE_IGNORE_FILES ${_file_name})

  _package_get_source_files(${_pkg} _srcs _pub_hdrs _priv_hdrs)
  list(APPEND CPACK_SOURCE_IGNORE_FILES ${_srcs} ${_pub_hdrs} ${_priv_hdrs})
endforeach()
list(APPEND CPACK_SOURCE_IGNORE_FILES "/.*build.*/;/CVS/;/\\\\.svn/;/\\\\.bzr/;/\\\\.hg/;/\\\\.hgignore;/\\\\.git/;\\\\.swp$;\\\\.#;/#;~")

include(CPack)
