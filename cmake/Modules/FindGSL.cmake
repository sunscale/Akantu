#===============================================================================
# @file   FindGSL.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Sun Oct 19 2014
#
# @brief  find package for gsl
#
# @section LICENSE
#
# Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory
# (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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

find_path(GSL_INCLUDE_PATH gsl_math.h
  PATHS ${GSL_DIR} ENV C_INCLUDE_PATH
  PATH_SUFFIXES gsl
  )

find_library(GSL_MAIN_LIBRARY NAME gsl
  PATHS ${GSL_DIR} ENV LIBRARY_PATH
  PATH_SUFFIXES lib
  )

find_library(GSL_BLAS_LIBRARY NAME gslcblas
  PATHS ${GSL_DIR} ENV LIBRARY_PATH
  PATH_SUFFIXES lib
)

mark_as_advanced(GSL_INCLUDE_PATH)
mark_as_advanced(GSL_MAIN_LIBRARY NAME)
mark_as_advanced(GSL_BLAS_LIBRARY NAME)

set(GSL_LIBRARIES ${GSL_MAIN_LIBRARY} ${GSL_BLAS_LIBRARY})


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GSL DEFAULT_MSG
  GSL_LIBRARIES GSL_INCLUDE_PATH)
