#===============================================================================
# @file   IOHelperCPack.cmake
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
#
# @date creation: Wed Oct 17 2012
# @date last modification: Mon Jul 28 2014
#
# @brief  configuration for CPack
#
# @section LICENSE
#
# Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# IOHelper is free  software: you can redistribute it and/or  modify it under the
# terms  of the  GNU Lesser  General Public  License as  published by  the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# IOHelper is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
# details.
#
# You should  have received  a copy  of the GNU  Lesser General  Public License
# along with IOHelper. If not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================

set(CPACK_GENERATOR "DEB;STGZ")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "guillaume.anciaux@epfl.ch")
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libc6, zlib1g")

if(CMAKE_SYSTEM_PROCESSOR MATCHES "i.86" OR CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64")
  if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64" CACHE STRING "Architecture of debian package generation")
  else()
    set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "i386" CACHE STRING "Architecture of debian package generation")
  endif()
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "ppc")
  set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "powerpc" CACHE STRING "Architecture of debian package generation")
else()
  set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "unknown" CACHE STRING "Architecture of debian package generation")
endif()

set(CPACK_DEBIAN_PACKAGE_DESCRIPTION "IOHelper library")
set(CPACK_PACKAGE_VENDOR "LSMS")
set(CPACK_PACKAGE_VERSION ${IOHelper_VERSION})

include(CPack)
