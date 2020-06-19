#===============================================================================
# @file   CMakeLists.txt
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Tue Sep 28 2010
# @date last modification: Mon Jan 18 2016
#
# @brief  configuration file for doxygen
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
find_package(Doxygen REQUIRED)

set(DOXYGEN_WARNINGS NO)
set(DOXYGEN_QUIET YES)
if(CMAKE_VERBOSE_MAKEFILE)
  set(DOXYGEN_WARNINGS YES)
  set(DOXYGEN_QUIET NO)
endif(CMAKE_VERBOSE_MAKEFILE)


set(DOXYGEN_INPUT ${CMAKE_CURRENT_BINARY_DIR}/akantu.dox)
set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/html)
set(DOXYGEN_OUTPUT ${DOXYGEN_OUTPUT_DIR}/index.html)

string(REGEX REPLACE ";" " " AKANTU_DOXYGEN_DEFINTIONS "${AKANTU_DEFINITIONS};DOXYGEN")
string(REGEX REPLACE ";" " " AKANTU_DOXYGEN_INCLUDE_DIRS "${AKANTU_INCLUDE_DIRS}")

make_directory(${DOXYGEN_OUTPUT_DIR})
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/akantu.dox.in
  ${CMAKE_CURRENT_BINARY_DIR}/akantu.dox)

add_custom_command(
  OUTPUT ${DOXYGEN_OUTPUT}
  COMMAND ${CMAKE_COMMAND} -E echo "${DOXYGEN_EXECUTABLE} ${DOXYGEN_INPUT}"
  COMMAND ${CMAKE_COMMAND} -E echo_append "Building akantu Documentation..."
  COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_INPUT}
  COMMAND ${CMAKE_COMMAND} -E echo "Done."
  DEPENDS ${DOXYGEN_INPUT}
  )

add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/latex/refman.pdf
  COMMAND ${CMAKE_COMMAND} -E echo "Building akantu RefMan..."
  COMMAND make -C ${CMAKE_CURRENT_BINARY_DIR}/latex
  COMMAND ${CMAKE_COMMAND} -E echo "Done."
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/latex/refman.tex
)

add_custom_target(doxygen-doc ALL DEPENDS ${DOXYGEN_OUTPUT})

add_custom_target(doxygen-doc-forced
  COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_INPUT}
  )

install(DIRECTORY ${DOXYGEN_OUTPUT_DIR}
  DESTINATION share/akantu-${AKANTU_VERSION}/doc)
