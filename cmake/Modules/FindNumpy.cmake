#===============================================================================
# @file   FindNumpy.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Mon Nov 16 2015
#
# @brief  The find_package file for numpy
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

find_package(PythonInterp)
#numpy includes
if(PYTHONINTERP_FOUND)
  set(_get_include "from __future__ import print_function; import numpy; print(numpy.get_include(), end='')")
  set(_get_version "from __future__ import print_function; import numpy; print(numpy.version.full_version, end=''),")

  execute_process(COMMAND
    ${PYTHON_EXECUTABLE} -c "${_get_include}"
    OUTPUT_VARIABLE _numpy_include_dir
    ERROR_QUIET
    RESULT_VARIABLE _res)

  if(_res EQUAL 0)
    set(NUMPY_INCLUDE_DIR "${_numpy_include_dir}" CACHE PATH "Include directory for numpy" FORCE)

    execute_process(COMMAND
      ${PYTHON_EXECUTABLE} -c "${_get_version}"
      OUTPUT_VARIABLE _numpy_version
      ERROR_QUIET
      RESULT_VARIABLE _res)
    if(_res EQUAL 0)
      set(NUMPY_VERSION "${_numpy_version}" CACHE STRING "Version of numpy")
    else()
      set(NUMPY_VERSION "NUMPY_VERSION-NOTFOUND" CACHE STRING "Version of numpy")
    endif()

    mark_as_advanced(NUMPY_INCLUDE_DIR NUMPY_VERSION)
  else()
    set(NUMPY_INCLUDE_DIR "NUMPY_INCLUDE_DIR-NOTFOUND" CACHE PATH "")
  endif()
endif()

#===============================================================================
include(FindPackageHandleStandardArgs)
if(CMAKE_VERSION VERSION_GREATER 2.8.12)
  find_package_handle_standard_args(Numpy
    REQUIRED_VARS NUMPY_INCLUDE_DIR
    VERSION_VAR NUMPY_VERSION)
else()
  find_package_handle_standard_args(Numpy DEFAULT_MSG
    NUMPY_INCLUDE_DIR NUMPY_VERSION)
endif()
