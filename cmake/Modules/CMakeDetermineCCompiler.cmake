#===============================================================================
# @file   CMakeDetermineCCompiler.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Sun Oct 19 2014
#
# @brief  CMake file to determine the compiler
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

macro(determine_compiler_version COMPILER)
  exec_program(${CMAKE_CXX_COMPILER}
    ARGS --version
    OUTPUT_VARIABLE _temp
    )

  set(${COMPILER}_COMPILER_VERSION "" CACHE STRING "Vesion of ${COMPILER} compiler.")
  string(REGEX MATCH "([0-9\\.]+)"
    ${COMPILER}_COMPILER_VERSION
    ${_temp}
    )

  mark_as_advanced(${COMPILER}_COMPILER_VERSION)
endmacro()

# Code from James Bigler (http://www.cmake.org/pipermail/cmake/2007-June/014460.html)
set(MANTA_COMPILER_NAME_REGEXPR "icc.*$")
if(NOT CMAKE_COMPILER_IS_GNUCC)
   # This regular expression also matches things like icc-9.1
   if(CMAKE_C_COMPILER MATCHES ${MANTA_COMPILER_NAME_REGEXPR})
     set(AKANTU_USING_ICC TRUE)
   endif(CMAKE_C_COMPILER MATCHES ${MANTA_COMPILER_NAME_REGEXPR})
else(NOT CMAKE_COMPILER_IS_GNUCC)
  set(AKANTU_USING_GNUCC TRUE)
endif(NOT CMAKE_COMPILER_IS_GNUCC)

set(MANTA_COMPILER_NAME_REGEXPR "icpc.*$")
if(NOT CMAKE_COMPILER_IS_GNUCXX)
  if(CMAKE_CXX_COMPILER MATCHES ${MANTA_COMPILER_NAME_REGEXPR})
    set(AKANTU_USING_ICPC TRUE)
    determine_compiler_version(INTEL)
    #else mvsc/clang/ibm/... ?
  endif(CMAKE_CXX_COMPILER MATCHES ${MANTA_COMPILER_NAME_REGEXPR})
else(NOT CMAKE_COMPILER_IS_GNUCXX)
  set(AKANTU_USING_GNUCXX TRUE)
  determine_compiler_version(GCC)
endif(NOT CMAKE_COMPILER_IS_GNUCXX)
