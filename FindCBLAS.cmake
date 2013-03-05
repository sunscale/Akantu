#===============================================================================
# @file   FindCBLAS.cmake
# @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
# @date   Thu Apr 19 16:44:00 2012
#
# @brief  The find_package file for the CBLAS library
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

find_library (CBLAS_LIBRARY cblas)

set(CBLAS_LIBRARIES ${CBLAS_LIBRARY} CACHE INTERNAL "Libraries for CBLAS" FORCE)

find_path (CBLAS_INCLUDE_DIR cblas.h)
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CBLAS DEFAULT_MSG CBLAS_INCLUDE_DIR CBLAS_LIBRARIES)