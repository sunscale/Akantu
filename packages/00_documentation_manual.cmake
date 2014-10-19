#===============================================================================
# @file   00_documentation_manual.cmake
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
#
# @date creation: Tue Jun 10 2014
# @date last modification: Thu Jul 03 2014
#
# @brief  Akantu's manual package
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

option(AKANTU_DOCUMENTATION_MANUAL  "Build manual." OFF)

set(AKANTU_DOCUMENTATION_MANUAL_DOCUMENTATION "
This package alows to compile the user manual in the build folder \\shellcode{build/doc/manual/manual.pdf}.

Under Ubuntu (14.04 LTS), the installation of the dependencies can be performed using the following command:
\\begin{command}
  > sudo apt-get install install rubber texlive texlive-science texlive-latex-extra
\\end{command}")