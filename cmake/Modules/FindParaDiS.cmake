#===============================================================================
# @file   FindParaDiS.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Sun Oct 19 2014
#
# @brief  find_package for paradis
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

find_library(PARADIS_LIBRARIES paradis
  PATHS ${PARADIS_DIR}
  PATH_SUFFIXES bin)

find_path(PARADIS_INCLUDE_PATH ParadisGen.h
  PATHS ${PARADIS_DIR}
  PATH_SUFFIXES include
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PARADIS DEFAULT_MSG
  PARADIS_LIBRARIES PARADIS_INCLUDE_PATH)

if (NOT PARADIS_FOUND)
  set(PARADIS_DIR "" CACHE PATH "Location of PARADIS library")
endif(NOT PARADIS_FOUND)