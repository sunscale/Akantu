#===============================================================================
# @file   core_cxx11.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Tue Feb 26 2013
# @date last modification: Wed Dec 16 2015
#
# @brief  C++11 addition to the core package
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

if(AKANTU_CXX11_FLAGS)
  package_declare(core_cxx11 ADVANCED
    DESCRIPTION "C++ 11 additions for Akantu core" DEFAULT ON
    COMPILE_FLAGS CXX "${AKANTU_CXX11_FLAGS}")

  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.6")
      set(AKANTU_CORE_CXX11 OFF CACHE BOOL "C++ 11 additions for Akantu core - not supported by the selected compiler" FORCE)
    endif()
  endif()
else()
  package_declare(core_cxx11 ADVANCED
    DESCRIPTION "C++ 11 additions for Akantu core"
    DEFAULT OFF
    NOT_OPTIONAL
    )
endif()

package_declare_sources(core_cxx11
  common/aka_point.hh
  common/aka_ball.cc
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


package_declare_documentation(core_cxx11
  "This option activates some features of the C++11 standard. This is usable with GCC>=4.7 or Intel>=13.")

