#===============================================================================
# @file   90_scotch.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Mon Nov 21 2011
# @date last modification: Thu Jul 10 2014
#
# @brief  package description for scotch
#
# @section LICENSE
#
# Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
package_declare(Scotch EXTERNAL
  DESCRIPTION "Add Scotch support in akantu"
  SYSTEM ON third-party/cmake/scotch.cmake)

package_declare_sources(Scotch
  mesh_utils/mesh_partition/mesh_partition_scotch.cc
  )

package_add_third_party_script_variable(Scotch
  SCOTCH_VERSION "5.1.12b")
package_add_third_party_script_variable(Scotch
  SCOTCH_ARCHIVE_HASH "MD5=e13b49be804755470b159d7052764dc0")
package_add_third_party_script_variable(Scotch
  SCOTCH_URL "https://gforge.inria.fr/frs/download.php/28978/scotch_${SCOTCH_VERSION}_esmumps.tar.gz")


package_declare_documentation(Scotch
  "This package enables the use the \\href{http://www.labri.fr/perso/pelegrin/scotch/}{Scotch}"
  "library in order to perform a graph partitioning leading to the domain"
  "decomposition used within \\akantu"
  ""
  "Under Ubuntu (14.04 LTS) the installation can be performed using the commands:"
  "\\begin{command}"
  "  > sudo apt-get install libscotch-dev"
  "\\end{command}"
  ""
  "If you activate the advanced option AKANTU\\_USE\\_THIRD\\_PARTY\\_SCOTCH"
  "the make system of akantu can automatically compile Scotch."
  ""
  "If the automated download fails due to a SSL access not supported by your"
  "version of CMake please download the file"
  "\\href{${SCOTCH_ARCHIVE}}{scotch\\_${SCOTCH_VERSION}\\_esmumps.tar.gz}"
  "and then place it in the directory \\shellcode{<akantu source>/third-party}"
 )

# if(SCOTCH_INCLUDE_DIR)
#   file(STRINGS ${SCOTCH_INCLUDE_DIR}/scotch.h SCOTCH_INCLUDE_CONTENT)
#   string(REGEX MATCH "_cplusplus" _match ${SCOTCH_INCLUDE_CONTENT})
#   if(_match)
#     set(AKANTU_SCOTCH_NO_EXTERN ON)
#     list(APPEND AKANTU_DEFINITIONS AKANTU_SCOTCH_NO_EXTERN)
#   else()
#     set(AKANTU_SCOTCH_NO_EXTERN OFF)
#   endif()
# endif()
