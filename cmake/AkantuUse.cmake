#===============================================================================
# @file   AkantuUse.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Tue Dec 07 2010
# @date last modification: Mon Aug 17 2015
#
# @brief  CMake file for the library
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
#===============================================================================
if (DEFINED CMAKE_PACKAGES_SYSTEM_LOADED)
  return()
endif()
set(CMAKE_PACKAGES_SYSTEM_LOADED TRUE)

function(package_is_activated pkg activated)
  string(TOUPPER ${pkg} _u_pkg)
  set(${activated} ${AKANTU_HAS_${_u_pkg}} PARENT_SCOPE)
endfunction()

function(package_get_include_dir pkg include_dir)
  string(TOUPPER ${pkg} _u_pkg)
  set(${include_dir} ${AKANTU_${_u_pkg}_INCLUDE_DIR} PARENT_SCOPE)
endfunction()

function(package_get_libraries pkg libs)
  string(TOUPPER ${pkg} _u_pkg)
  set(${libs} ${AKANTU_${_u_pkg}_LIBRARIES} PARENT_SCOPE)
endfunction()

function(package_get_compile_flags pkg lang flags)
  string(TOUPPER ${pkg} _u_pkg)
  set(${flags} ${AKANTU_${_u_pkg}_COMPILE_${lang}_FLAGS} PARENT_SCOPE)
endfunction()


function(get_target_list_of_associated_files tgt files)
  get_target_property(_type ${tgt} TYPE)
  if(_type STREQUAL "SHARED_LIBRARY"
      OR _type STREQUAL "STATIC_LIBRARY"
      OR _type STREQUAL "MODULE_LIBRARY"
      OR _type STREQUAL "EXECUTABLE")
    get_target_property(_srcs ${tgt} SOURCES)
    set(_dep_ressources)
    foreach(_file ${_srcs})
      list(APPEND _dep_ressources ${CMAKE_CURRENT_SOURCE_DIR}/${_file})
    endforeach()
  else()
    get_target_property(_dep_ressources ${tgt} RESSOURCES)
  endif()

  set(${files} ${_dep_ressources} PARENT_SCOPE)
endfunction()
