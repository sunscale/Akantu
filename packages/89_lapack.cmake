#===============================================================================
# @file   89_lapack.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Fri Oct 19 2012
# @date last modification: Thu Jun 12 2014
#
# @brief  package description for lapack support
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

add_optional_external_package(LAPACK "Use LAPACK for arithmetic operations" OFF LANGUAGE Fortran)
mark_as_advanced(AKANTU_USE_LAPACK)

set(AKANTU_LAPACK_DOCUMENTATION "
This package provides access to a LAPACK implementation.

Under Ubuntu (14.04 LTS), the installation can be performed using the following command:
\\begin{command}
  > sudo apt-get install libatlas-base-dev
\\end{command}
" )
