#===============================================================================
# @file   scalapack.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Fri Oct 19 2012
# @date last modification: Mon Mar 30 2015
#
# @brief  package description for mumps support
#
# @section LICENSE
#
# Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
# Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
# Solides)
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

package_declare(ScaLAPACK EXTERNAL
  DESCRIPTION "Add ScaLAPACK support in akantu"
  SYSTEM OFF third-party/cmake/scalapack.cmake
  DEPENDS MPI
  )


package_add_third_party_script_variable(ScaLAPACK
  SCALAPACK_VERSION "2.0.2")
package_add_third_party_script_variable(ScaLAPACK
  SCALAPACK_ARCHIVE "http://www.netlib.org/scalapack/scalapack-${SCALAPACK_VERSION}.tgz")
package_add_third_party_script_variable(ScaLAPACK
  SCALAPACK_ARCHIVE_HASH_2.0.2 "MD5=2f75e600a2ba155ed9ce974a1c4b536f")
