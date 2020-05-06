#===============================================================================
# @file   documentation_doxygen.cmake
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Wed Jun 11 2014
# @date last modification: Mon Jan 18 2016
#
# @brief  Doxygen documentation of the code
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
package_declare(documentation_developer_manual
  DESCRIPTION "Build source documentation using Sphinx/Doxygen.")

package_declare_documentation(documentation_developer_manual
  "This generates the Doxygen documantation of the source code."
  "It depends on:"
  "\\begin{itemize}"
  "\\item \\href{http://www.stack.nl/~dimitri/doxygen/}{Doxygen} an automated source code documentations system."
  "\\item Optional: \\href{http://www.graphviz.org/}{Graphviz} to generate the dependencies graph"
  "\\end{itemize}"
  ""
  "Under Ubuntu (14.04 LTS), the installation of the dependencies can be performed using the following command:"
  "\\begin{command}"
  "  > sudo apt-get install doxygen"
  "  > sudo apt-get install graphviz"
  "\\end{command}"
  )

package_set_package_system_dependency(documentation_developer_manual deb-src
  python3-sphinx python3-breathe doxygen graphviz)

package_declare_extra_files_to_package(documentation_developer_manual
  PROJECT doc/dev-doc/doxygen/akantu.dox.in
          doc/dev-doc/sphinx/conf.py.in
          doc/dev-doc/sphinx/index.rst
          cmake/Modules/FindSphinx.cmake
  )
