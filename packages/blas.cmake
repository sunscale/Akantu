#===============================================================================
# @file   90_blas.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Fri Oct 19 2012
# @date last modification: Tue Jun 24 2014
#
# @brief  package description for blas support
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

package_declare(BLAS EXTERNAL
  DESCRIPTION "Use BLAS for arithmetic operations"
  EXTRA_PACKAGE_OPTIONS LANGUAGE Fortran)

set(AKANTU_USE_BLAS_VENDOR "Generic" CACHE STRING "Version of blas to use")
mark_as_advanced(AKANTU_USE_BLAS_VENDOR)
set_property(CACHE AKANTU_USE_BLAS_VENDOR PROPERTY STRINGS
  Goto
  ATLAS
  PhiPACK
  CXML
  DXML
  SunPerf
  SCSL
  SGIMATH
  IBMESSL
  Intel10_32
  Intel10_64lp
  Intel10_64lp_seq
  Intel
  ACML
  ACML_MP
  ACML_GPU
  Apple
  NAS
  Generic
  )

set(ENV{BLA_VENDOR} ${AKANTU_USE_BLAS_VENDOR})

if(BLAS_mkl_core_LIBRARY)
  set(AKANTU_USE_BLAS_MKL CACHE INTERNAL "" FORCE)
endif()

package_declare_documentation(BLAS
  "This package provides access to a BLAS implementation."
  ""
  "Under Ubuntu (14.04 LTS), the installation can be performed using the following command:"
  "\\begin{command}"
  "  > sudo apt-get install libatlas-base-dev"
  "\\end{command}"
  )