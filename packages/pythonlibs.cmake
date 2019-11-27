#===============================================================================
# @file   pythonlibs.cmake
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Fri Sep 03 2010
# @date last modification: Fri Jan 22 2016
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
package_declare(PythonLibsNew EXTERNAL DESCRIPTION "Akantu's python interface"
  DEPENDS numpy
  EXTRA_PACKAGE_OPTIONS ARGS ${AKANTU_PREFERRED_PYTHON_VERSION} PREFIX PYTHON FOUND PYTHONLIBS_FOUND
  )

package_set_package_system_dependency(PythonLibs deb libpython3)
package_set_package_system_dependency(PythonLibs deb-src libpython3-dev)

package_declare_documentation(PythonLibs
  "This package is a dependency of the python interface"
  ""
  "Under Ubuntu (14.04 LTS) the installation can be performed using the commands:"
  "\\begin{command}"
  "  > sudo apt-get install libpython3-dev"
  "\\end{command}"
  ""
)
