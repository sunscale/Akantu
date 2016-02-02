#===============================================================================
# @file   FindEPSN.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Sun Oct 19 2014
#
# @brief  The find_package file for EPSN
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

#===============================================================================
find_path(EPSN_DIR EPSNConfig.cmake
  PATHS $ENV{EPSN_TOP}
  )


if(EPSN_DIR)
  include(${EPSN_DIR}/EPSNConfig.cmake)
  set(EPSN_LIB_PATH ${EPSN_DIR}/lib
    ${EPSN_LIBRARIES_DIR}
    )
  find_library(EPSN_COMMON_LIBRARY epsn_common
    PATHS ${EPSN_LIB_PATH}
    )
  find_library(EPSN_SIMULATION_LIBRARY epsn_simulation
    PATHS ${EPSN_LIB_PATH}
    )
  include(${EPSN_DIR}/EPSNLibraryDepends.cmake)

  set(EPSN_LIBRARIES ${EPSN_SIMULATION_LIBRARY} ${EPSN_COMMON_LIBRARY})
endif(EPSN_DIR)

#===============================================================================
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(EPSN DEFAULT_MSG
  EPSN_LIBRARIES EPSN_INCLUDE_DIR)
