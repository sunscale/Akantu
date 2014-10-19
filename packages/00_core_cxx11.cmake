#===============================================================================
# @file   00_core_cxx11.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Mon May 06 2013
# @date last modification: Thu Jul 03 2014
#
# @brief  C++11 addition to the core package
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

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag (-std=c++0x HAVE_NEW_STD)

set(AKANTU_CORE_CXX11_FILES
  common/aka_point.hh
  common/aka_ball.cc
  common/aka_ci_string.hh
  common/aka_plane.hh
  common/aka_polytope.hh
  common/aka_ball.hh
  common/aka_timer.hh
  common/aka_tree.hh
  common/aka_bounding_box.hh
  common/aka_bounding_box.cc
  common/aka_geometry.hh
  common/aka_geometry.cc
  model/solid_mechanics/solid_mechanics_model_element.hh
  )


if(HAVE_NEW_STD)
  option(AKANTU_CORE_CXX11 "C++ 11 additions for Akantu core" ON)

  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.6")
      set(AKANTU_CORE_CXX11 OFF CACHE BOOL "C++ 11 additions for Akantu core - not supported by the selected compiler" FORCE)
    elseif(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "4.8")
      if(AKANTU_CORE_CXX11)
	add_flags(cxx -DBOOST_RESULT_OF_USE_TR1)
      else()
	remove_flags(cxx -DBOOST_RESULT_OF_USE_TR1)
      endif()
    endif()
  endif()

  if(AKANTU_CORE_CXX11)
    add_flags(cxx "-std=c++0x")
  else()
    remove_flags(cxx "-std=c++0x")
  endif()
else()
  set(AKANTU_CORE_CXX11 OFF CACHE BOOL "core package for Akantu" FORCE)
  remove_flags(cxx "-std=c++0x")
endif()

mark_as_advanced(AKANTU_CORE_CXX11)

set(AKANTU_CORE_CXX11_DOCUMENTATION "
This option activates some features of the C++11 standard. This is usable with GCC>=4.7 or Intel>=13.")

