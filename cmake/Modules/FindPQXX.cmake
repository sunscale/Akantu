#===============================================================================
# @file   FindPQXX.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Sun Oct 19 2014
#
# @brief  The find_package file for PostgreSQL C++ library
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

find_package(PostgreSQL REQUIRED)
if(POSTGRESQL_FOUND)
  find_library(PQXX_LIBRARY NAMES pqxx
    HINTS ${PQXX_DIR} ENV PQXX_DIR
    DOC "Location of libpqxx library"
    )

  find_path(PQXX_HEADER_DIR NAMES pqxx/pqxx
    HINTS ${PQXX_DIR} ENV PQXX_DIR
    DOC "Path to pqxx/pqxx header file. Do not include the 'pqxx' directory in this value."
    )

  set(PQXX_INCLUDE_DIR "${PQXX_HEADER_DIR};${POSTGRESQL_INCLUDE_DIR}" CACHE STRING "Include directories for PostgreSQL C++ library"  FORCE)
  set(PQXX_LIBRARIES "${PQXX_LIBRARY};${POSTGRESQL_LIBRARIES}" CACHE STRING "Link libraries for PostgreSQL C++ interface" FORCE)

  mark_as_advanced(PQXX_HEADER_DIR)
  mark_as_advanced(PQXX_INCLUDE_DIR)
  mark_as_advanced(PQXX_LIBRARY)
  mark_as_advanced(PQXX_LIBRARIES)
endif()

#===============================================================================
if(NOT PQXX_FOUND)
  set(PQXX_DIR "" CACHE PATH "Help to find the location of pqxx library.")
  mark_as_advanced(PQXX_FOUND)
endif()

#===============================================================================
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PQXX DEFAULT_MSG PQXX_LIBRARY PQXX_HEADER_DIR)

