#===============================================================================
# @file   FindCppArray.cmake
# @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
# @date   Thu Apr 19 16:48:00 2012
#
# @brief  The find_package file for cpp-array library
#
# @section LICENSE
#
# Copyright (©) 2010-2012 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

find_path (CPPARRAY_INCLUDE_DIR expr.hpp PATH_SUFFIXES array)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CppArray DEFAULT_MSG CPPARRAY_INCLUDE_DIR)

