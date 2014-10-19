#===============================================================================
# @file   80_nlopt.cmake
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Thu Jun 05 2014
# @date last modification: Thu Sep 18 2014
#
# @brief  package for the opitmization library NLopt
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

option(AKANTU_USE_THIRD_PARTY_NLOPT "Automatic download of the NLOPT library" ON)
option(AKANTU_USE_NLOPT "Use NLOPT library")
mark_as_advanced(AKANTU_USE_THIRD_PARTY_NLOPT AKANTU_USE_NLOPT)

set(NLOPT_VERSION "2.4.2")
set(NLOPT_ARCHIVE "${PROJECT_SOURCE_DIR}/third-party/nlopt-${NLOPT_VERSION}.tar.gz")
set(NLOPT_ARCHIVE_HASH "MD5=d0b8f139a4acf29b76dbae69ade8ac54")
if(NOT EXISTS ${NLOPT_ARCHIVE})
  set(NLOPT_ARCHIVE "http://ab-initio.mit.edu/nlopt/nlopt-${NLOPT_VERSION}.tar.gz")
endif()

if (AKANTU_USE_THIRD_PARTY_NLOPT AND AKANTU_USE_NLOPT)
  set(NLOPT_CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR> --enable-shared --with-cxx)
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

  set(NLOPT_LIBRARIES   ${NLOPT_DIR}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}nlopt_cxx${CMAKE_SHARED_LIBRARY_SUFFIX} CACHE PATH "Libraries for NLopt" FORCE)
  set(NLOPT_INCLUDE_DIR ${NLOPT_DIR}/include CACHE PATH "Include directories for NLopt" FORCE)

  mark_as_advanced(NLOPT_INCLUDE_DIR NLOPT_LIBRARIES)

  list(APPEND AKANTU_EXTERNAL_LIBRARIES ${NLOPT_LIBRARIES})
  list(APPEND AKANTU_EXTERNAL_LIB_INCLUDE_DIR ${NLOPT_INCLUDE_DIR})
  set(AKANTU_NLOPT_INCLUDE_DIR ${NLOPT_INCLUDE_DIR})
  set(AKANTU_NLOPT_LIBRARIES ${NLOPT_LIBRARIES})
  list(APPEND AKANTU_OPTION_LIST NLOPT)
  set(NLOPT_FOUND TRUE CACHE INTERNAL "" FORCE)
  set(AKANTU_NLOPT ON)

  list(APPEND AKANTU_EXTRA_TARGET_DEPENDENCIES NLopt)
else()
  add_optional_external_package(NLopt "Use NLOPT library" OFF)
endif()


set(AKANTU_NLOPT_DOCUMENTATION "
This package enable the use of the optimization alogorithm library \\href{http://ab-initio.mit.edu/wiki/index.php/NLopt}{NLopt}.
Since there are no packaging for the common operating system by default \\akantu compiles it as a third-party project. This behavior can be modified with the option \\code{AKANTU\\_USE\\_THIRD\\_PARTY\\_NLOPT}.

If the automated download fails please download \\href{http://ab-initio.mit.edu/nlopt/nlopt-${NLOPT_VERSION}.tar.gz}{nlopt-${NLOPT_VERSION}.tar.gz} and place it in \\shellcode{<akantu source>/third-party} download.
")
