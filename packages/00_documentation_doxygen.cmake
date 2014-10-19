#===============================================================================
# @file   00_documentation_doxygen.cmake
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
#
# @date creation: Tue Jun 10 2014
# @date last modification: Tue Jun 24 2014
#
# @brief  Doxygen documentation of the code
#
# @section LICENSE
#
# Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

option(AKANTU_DOCUMENTATION_DOXYGEN "Build source documentation using Doxygen." OFF)

set(AKANTU_DOCUMENTATION_DOXYGEN_DOCUMENTATION
"
This generates the Doxygen documantation of the source code.
It depends on:
\\begin{itemize}
\\item \\href{http://www.stack.nl/~dimitri/doxygen/}{Doxygen} an automated source code documentations system.
\\item Optional: \\href{http://www.graphviz.org/}{Graphviz} to generate the dependencies graph
\\end{itemize}

Under Ubuntu (14.04 LTS), the installation of the dependencies can be performed using the following command:
\\begin{command}
  > sudo apt-get install doxygen
  > sudo apt-get install graphviz
\\end{command}
")