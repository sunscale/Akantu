#===============================================================================
# @file   numpy.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Tue Nov 29 2011
# @date last modification: Tue Jan 19 2016
#
# @brief  package description for the python library
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

package_declare(Numpy EXTERNAL DESCRIPTION "Akantu's numpy dependance check")

package_declare_documentation(Numpy
  "This package allows to wrap Akantu arrays to numpy arrays"
  ""
  "Under Ubuntu (14.04 LTS) the installation can be performed using the commands:"
  "\\begin{command}"
  "  > sudo apt-get install python-numpy"
  "\\end{command}"
  ""
)

package_set_package_system_dependency(Numpy deb python-numpy)
package_set_package_system_dependency(Numpy deb-src python-numpy)

package_declare_extra_files_to_package(Numpy
  PROJECT cmake/Modules/FindNumpy.cmake)
