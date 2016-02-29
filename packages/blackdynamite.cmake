#===============================================================================
# @file   blackdynamite.cmake
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Fri Mar 15 2013
# @date last modification: Mon Jan 18 2016
#
# @brief  package description for BlackDynamite support
#
# @section LICENSE
#
# Copyright  (©)  2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de Lausanne)
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
  SYSTEM OFF third-party/cmake/blackdynamite.cmake
  EXTRA_PACKAGE_OPTIONS FOUND BlackDynamite_FOUND)

set(_version master)

package_add_third_party_script_variable(BlackDynamite
  BLACKDYNAMITE_VERSION "${_version}")
package_add_third_party_script_variable(BlackDynamite
  BLACKDYNAMITE_GIT "git@lsmssrv1.epfl.ch:blackdynamite.git")
package_add_third_party_script_variable(BlackDynamite
  BLACKDYNAMITE_ARCHIVE "blackdynamite-${_version}.tar.gz")

package_declare_extra_files_to_package(BlackDynamite
  PROJECT third-party/cmake/blackdynamite.cmake
  )
