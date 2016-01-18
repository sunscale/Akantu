#===============================================================================
# @file   nlopt.cmake
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Thu Jun 05 2014
# @date last modification: Mon Mar 30 2015
#
# @brief  package for the opitmization library NLopt
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

package_declare(nlopt EXTERNAL
  DESCRIPTION "Use NLOPT library"
  SYSTEM OFF)

package_use_system(nlopt _use_system)

if(NOT ${_use_system})
  package_get_option_name(nlopt _option_name)

  if(${_option_name})
    set(NLOPT_VERSION "2.4.2")
    set(NLOPT_ARCHIVE "${PROJECT_SOURCE_DIR}/third-party/nlopt-${NLOPT_VERSION}.tar.gz")
    set(NLOPT_ARCHIVE_HASH "MD5=d0b8f139a4acf29b76dbae69ade8ac54")
    if(NOT EXISTS ${NLOPT_ARCHIVE})
      set(NLOPT_ARCHIVE "http://ab-initio.mit.edu/nlopt/nlopt-${NLOPT_VERSION}.tar.gz")
    endif()

    string(TOUPPER "${CMAKE_BUILD_TYPE}" _u_build_type)
    set(NLOPT_CONFIGURE_COMMAND CXX=${CMAKE_CXX_COMPILER} CXX_FLAGS=${CMAKE_CXX_FLAGS_${_u_build_type}} <SOURCE_DIR>/configure
      --prefix=<INSTALL_DIR> --enable-shared --with-cxx --without-threadlocal
      --without-guile --without-python --without-octave --without-matlab)
  set(NLOPT_DIR ${PROJECT_BINARY_DIR}/third-party)

    include(ExternalProject)

    ExternalProject_Add(NLopt
      PREFIX ${NLOPT_DIR}
      URL ${NLOPT_ARCHIVE}
      URL_HASH ${NLOPT_ARCHIVE_HASH}
      CONFIGURE_COMMAND ${NLOPT_CONFIGURE_COMMAND}
      BUILD_COMMAND make
      INSTALL_COMMAND make install
      )

    set_third_party_shared_libirary_name(NLOPT_LIBRARIES nlopt_cxx)
    set(NLOPT_INCLUDE_DIR ${NLOPT_DIR}/include CACHE PATH "Include directories for NLopt" FORCE)
    mark_as_advanced(NLOPT_INCLUDE_DIR)

    package_set_libraries(nlopt ${NLOPT_LIBRARIES})
    package_set_include_dir(nlopt ${NLOPT_INCLUDE_DIR})

    package_add_extra_dependency(nlopt NLopt)
  endif()
else()
  package_rm_extra_dependency(nlopt NLopt)
endif()

package_declare_documentation(nlopt
  "This package enable the use of the optimization alogorithm library \\href{http://ab-initio.mit.edu/wiki/index.php/NLopt}{NLopt}."
  "Since there are no packaging for the common operating system by default \\akantu compiles it as a third-party project. This behavior can be modified with the option \\code{AKANTU\\_USE\\_THIRD\\_PARTY\\_NLOPT}."
  ""
  "If the automated download fails please download \\href{http://ab-initio.mit.edu/nlopt/nlopt-${NLOPT_VERSION}.tar.gz}{nlopt-${NLOPT_VERSION}.tar.gz} and place it in \\shellcode{<akantu source>/third-party} download."
  )

package_declare_extra_files_to_package(nlopt
  PROJECT cmake/Modules/FindNLopt.cmake)
