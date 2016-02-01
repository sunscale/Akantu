#===============================================================================
# @file   FindPostgreSQL.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Sun Oct 19 2014
#
# @brief  The find_package file for PostgreSQL C library
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

find_library(POSTGRESQL_LIBRARIES NAMES pq
  PATH_SUFFIXES pgsql postgresql
  PATH ENV POSTGRESQL_DIR
  )

find_path(POSTGRESQL_INCLUDE_DIR NAMES libpq-fe.h
  PATH_SUFFIXES pgsql postgresql
  PATHS ENV POSTGRESQL_DIR
  )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(POSTGRESQL DEFAULT_MSG
  POSTGRESQL_LIBRARIES POSTGRESQL_INCLUDE_DIR)

mark_as_advanced(POSTGRESQL_INCLUDE_DIR)
mark_as_advanced(POSTGRESQL_LIBRARIES)
