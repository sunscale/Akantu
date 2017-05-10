#===============================================================================
# @file   blas.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Tue Oct 16 2012
# @date last modification: Mon Jan 18 2016
#
# @brief  package description for blas support
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

package_declare(BLAS EXTERNAL
  DESCRIPTION "Use BLAS for arithmetic operations"
  EXTRA_PACKAGE_OPTIONS LANGUAGE Fortran
  SYSTEM ON third-party/cmake/blas.cmake)

package_add_third_party_script_variable(BLAS BLAS_ARCHIVE "http://www.netlib.org/blas/blas-3.5.0.tgz")
package_add_third_party_script_variable(BLAS BLAS_VERSION "3.5.0")

set(_default_blas $ENV{BLA_VENDOR})
if(NOT _default_blas)
  set(_default_blas Generic)
endif()

set(AKANTU_USE_BLAS_VENDOR "${_default_blas}" CACHE STRING "Version of blas to use")
mark_as_advanced(AKANTU_USE_BLAS_VENDOR)
set_property(CACHE AKANTU_USE_BLAS_VENDOR PROPERTY STRINGS
  ACML
  ACML_GPU
  ACML_MP
  ATLAS
  Apple
  CXML
  DXML
  Generic
  Goto
  IBMESSL
  Intel
  Intel10_32
  Intel10_64lp
  Intel10_64lp_seq
  NAS
  OpenBLAS
  PhiPACK
  SCSL
  SGIMATH
  SunPerf
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

package_set_package_system_dependency(BLAS deb libblas3)
package_set_package_system_dependency(BLAS deb-src libblas3)

package_declare_extra_files_to_package(BLAS
  PROJECT
    third-party/cmake/blas.cmake
    third-party/blas_3.5.0_make.inc.cmake
  )
