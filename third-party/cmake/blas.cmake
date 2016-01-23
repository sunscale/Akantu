#===============================================================================
# @file   blas.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Wed Nov 11 2015
# @date last modification: Wed Nov 11 2015
#
# @brief  build script for netlib-blas
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

set(BLAS_DIR ${PROJECT_BINARY_DIR}/third-party)

configure_file(
  ${PROJECT_SOURCE_DIR}/third-party/blas_${BLAS_VERSION}_make.inc.cmake
  ${BLAS_DIR}/blas_make.inc @ONLY)

file(MAKE_DIRECTORY ${BLAS_DIR}/lib)

ExternalProject_Add(netlib-blas
  PREFIX ${BLAS_DIR}
  URL ${BLAS_ARCHIVE}
  CONFIGURE_COMMAND cmake -E copy ${BLAS_DIR}/blas_make.inc make.inc
  BUILD_IN_SOURCE 1
  BUILD_COMMAND ${CMAKE_MAKE_PROGRAM}
  INSTALL_COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SHARED_LIBRARY_PREFIX}blas${CMAKE_SHARED_LIBRARY_SUFFIX} <INSTALL_DIR>/lib
  LOG_DOWNLOAD 1
  LOG_CONFIGURE 1
  LOG_BUILD 1
  LOG_INSTALL 1
  )

package_add_extra_dependency(BLAS netlib-blas)

set_third_party_shared_libirary_name(BLAS_LIBRARIES blas)
