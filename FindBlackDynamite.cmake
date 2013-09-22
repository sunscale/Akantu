#===============================================================================
# @file   FindIOHelper.cmake
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @date   Mon Jun 27 16:29:57 2010
#
# @brief  The find_package file for libQVIEW
#
# @section LICENSE
#
# Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

find_library(BlackDynamite_LIBRARIES NAME blackdynamite
  PATHS ${BlackDynamite_DIR} 
  PATH_SUFFIXES lib
  )
#===============================================================================
find_path(BlackDynamite_INCLUDE_PATH pusher.hh
  PATHS ${BlackDynamite_DIR} ENV C_INCLUDE_PATH
  PATH_SUFFIXES include src
  )
#===============================================================================
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(BlackDynamite DEFAULT_MSG
  BlackDynamite_LIBRARIES BlackDynamite_INCLUDE_PATH)

#===============================================================================
if(NOT BlackDynamite_FOUND)
  set(BlackDynamite_DIR "" CACHE PATH "Location of BlackDynamite library.")
endif(NOT BlackDynamite_FOUND)

