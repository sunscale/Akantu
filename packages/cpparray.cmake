#===============================================================================
# @file   cpparray.cmake
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Fri Mar 15 2013
# @date last modification: Mon Mar 30 2015
#
# @brief  package description for cpp_array project
#
# @section LICENSE
#
# Copyright  (©)  2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de Lausanne)
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

package_declare(CppArray EXTERNAL
  DESCRIPTION "Use cpp-array library"
  SYSTEM OFF)

package_use_system(CppArray _use_system)
package_get_option_name(CppArray _option_name)

if(NOT ${_use_system})
  if(${_option_name})
    if(TARGET cpparray)
      return()
    endif()

    set(CPPARRAY_DIR ${PROJECT_BINARY_DIR}/third-party)

    include(ExternalProject)
    ExternalProject_Add(cpparray
      PREFIX ${CPPARRAY_DIR}
      GIT_REPOSITORY https://code.google.com/p/cpp-array.git
      CMAKE_ARGS <SOURCE_DIR>/cpp-array
      CMAKE_CACHE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR> -Dcpp-array_TESTS:BOOL=OFF -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE} -DCMAKE_CXX_COMPILER:PATH=${CMAKE_CXX_COMPILER} -DCMAKE_Fortran_COMPILER:PATH=${CMAKE_Fortran_COMPILER}
      )

    set(CPPARRAY_INCLUDE_DIR ${CPPARRAY_DIR}/include CACHE PATH "" FORCE)

    package_set_include_dir(CppArray ${CPPARRAY_INCLUDE_DIR})
    package_add_extra_dependency(CppArray cpparray)
  endif()
endif()

package_declare_documentation(CppArray
  "This package provides access to the \\href{https://code.google.com/p/cpp-array/}{cpp-array}"
  "open-source project. If internet is accessible when configuring the project (during cmake call)"
  "this package will be auto-downloaded."
  )
