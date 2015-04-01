#===============================================================================
# @file   blackdynamite.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date   Tue Nov 29 15:16:35 2011
#
# @brief  package description for BlackDynamite support
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
package_declare(BlackDynamite EXTERNAL
  DESCRIPTION "Use BlackDynamite library"
  SYSTEM OFF)

package_use_system(BlackDynamite _use_system)

if(NOT ${_use_system})
  package_get_option_name(BlackDynamite _option_name)

  if(${_option_name})
    set(BLACKDYNAMITE_URL "svn+ssh://lsmssrv1.epfl.ch/space/repositories/SimulPack/BlackDynamite")

    include(ExternalProject)

    ExternalProject_Add(blackdynamite
      PREFIX ${PROJECT_BINARY_DIR}/third-party
      SVN_REPOSITORY ${BLACKDYNAMITE_URL}
      CMAKE_ARGS <SOURCE_DIR>/
      CMAKE_CACHE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR> -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE} -DCMAKE_CXX_COMPILER:PATH=${CMAKE_CXX_COMPILER}
      BUILD_COMMAND make
      INSTALL_COMMAND make install
      )

    set_third_party_shared_libirary_name(BLACKDYNAMITE_LIBRARIES blackdynamite)
    package_set_libraries(BlackDynamite ${BLACKDYNAMITE_LIBRARIES})
    package_set_include_dir(BlackDynamite ${PROJECT_BINARY_DIR}/third-party/include/blackdynamite)

    package_add_extra_dependency(BlackDynamite blackdynamite)
  endif()
endif()
