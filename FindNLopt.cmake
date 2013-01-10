#===============================================================================
# @file   FindNLopt.cmake
# @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
# @date   Thu Apr 19 16:40:00 2012
#
# @brief  The find_package file for NLopt optimization library
#
# @section LICENSE
#
# Copyright (©) 2010-2012 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#===============================================================================
find_library(NLOPT_LIBRARIES NAMES nlopt_cxx
  PATHS ${NLOPT_DIR}
  PATH_SUFFIXES lib
  )

find_path(NLOPT_INCLUDE_DIR nlopt.hpp
  PATHS ${NLOPT_DIR}
  PATH_SUFFIXES include
  )



#===============================================================================
mark_as_advanced(NLOPT_LIBRARIES)
mark_as_advanced(NLOPT_INCLUDE_DIR)

#===============================================================================
if(NOT NLOPT_FOUND)
  set(NLOPT_DIR "" CACHE PATH "Location of NLOPT source directory.")
endif()

#===============================================================================
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NLopt DEFAULT_MSG NLOPT_LIBRARIES NLOPT_INCLUDE_DIR)

