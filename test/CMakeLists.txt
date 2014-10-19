#===============================================================================
# @file   CMakeLists.txt
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Fri Sep 03 2010
# @date last modification: Thu Jul 03 2014
#
# @brief  configuration for tests
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
# @section DESCRIPTION
#
#===============================================================================

include_directories(
  ${AKANTU_INCLUDE_DIRS}
  ${AKANTU_EXTERNAL_LIB_INCLUDE_DIR}
  )

set(AKANTU_TESTS_FILES CACHE INTERNAL "")
#===============================================================================
# List of tests
#===============================================================================
add_akantu_test(test_common "Test the common part of Akantu")
add_akantu_test(test_static_memory "Test static memory")
add_akantu_test(test_fem "Test finite element functionalties")
add_akantu_test(test_mesh_utils "Test mesh utils")
add_akantu_test(test_model "Test model objects")
add_akantu_test(test_solver "Test solver function" PACKAGE mumps)
add_akantu_test(test_io "Test the IO modules")
add_akantu_test(test_contact "Test the contact part of Akantu" PACKAGE contact) 

add_akantu_test(test_surface_extraction "Test mesh utils surface extraction")
add_akantu_test(test_synchronizer "Test synchronizers" PACKAGE parallel)

file(GLOB_RECURSE __all_files "*.*")
set(_all_files)
foreach(_file ${__all_files})
  if("${_file}" MATCHES  "${PROJECT_SOURCE_DIR}/test")
    file(RELATIVE_PATH __file ${PROJECT_SOURCE_DIR} ${_file})
    list(APPEND _all_files ${__file})
  endif()
endforeach()

file(GLOB_RECURSE __all_files "*CMakeLists.txt")
foreach(_file ${__all_files})
  if("${_file}" MATCHES  "${PROJECT_SOURCE_DIR}/test")
    file(RELATIVE_PATH __file ${PROJECT_SOURCE_DIR} ${_file})
    list(APPEND _cmakes ${__file})
  endif()
endforeach()
list(APPEND AKANTU_TESTS_FILES ${_cmakes})

foreach(_file ${_all_files})
  list(FIND AKANTU_TESTS_FILES ${_file} _ret)
  if(_ret EQUAL -1)
    list(APPEND AKANTU_TESTS_EXCLUDE_FILES /${_file})
  endif()
endforeach()

set(AKANTU_TESTS_EXCLUDE_FILES ${AKANTU_TESTS_EXCLUDE_FILES} PARENT_SCOPE)

#foreach(f ${AKANTU_TESTS_EXCLUDE_FILES})
#  message(${f})
#endforeach()
