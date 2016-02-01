#===============================================================================
# @file   lapack.cmake
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Fri Oct 19 2012
# @date last modification: Mon Jan 18 2016
#
# @brief  package description for lapack support
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

package_declare(LAPACK EXTERNAL
  DESCRIPTION "Use LAPACK for arithmetic operations"
  EXTRA_PACKAGE_OPTIONS LANGUAGE Fortran)

package_declare_documentation(LAPACK
  "This package provides access to a LAPACK implementation."
  ""
  "Under Ubuntu (14.04 LTS), the installation can be performed using the following command:"
  "\\begin{command}"
  "  > sudo apt-get install libatlas-base-dev"
  "\\end{command}"
  )

package_set_package_system_dependency(LAPACK deb liblapack3)
package_set_package_system_dependency(LAPACK deb-src liblapack-dev)
