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

package_on_enabled_script(PythonLibsNew
  "set(PYTHON_MODULE_PREFIX \${PYTHON_MODULE_PREFIX} CACHE INTERNAL \"\")
set(PYTHON_MODULE_EXTENSION \${PYTHON_MODULE_EXTENSION} CACHE INTERNAL \"\")
set(PYTHON_VERSION_MAJOR \${PYTHON_VERSION_MAJOR} CACHE INTERNAL \"\")
set(PYTHON_VERSION_MINOR \${PYTHON_VERSION_MINOR} CACHE INTERNAL \"\")
set(PYTHON_SITE_PACKAGES \${PYTHON_SITE_PACKAGES} CACHE INTERNAL \"\")
")


if(AKANTU_PREFERRED_PYTHON_VERSION VERSION_GREATER 2.9)
  package_set_package_system_dependency(PythonLibsNew deb libpython3)
  package_set_package_system_dependency(PythonLibsNew deb-src libpython3-dev)
endif()

package_declare_documentation(PythonLibsNew
  "This package is a dependency of the python interface"
  ""
  "Under Ubuntu (14.04 LTS) the installation can be performed using the commands:"
  "\\begin{command}"
  "  > sudo apt-get install libpython3-dev"
  "\\end{command}"
  ""
)
