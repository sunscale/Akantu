#===============================================================================
# @file   boost.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Fri Sep 03 2010
# @date last modification: Fri Dec 11 2015
#
# @brief  package handling the dependencies to boost
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

set(Boost_NO_BOOST_CMAKE ON CACHE BOOL "" FORCE)

package_declare(Boost EXTERNAL
  NOT_OPTIONAL
  DESCRIPTION "Package handling boost components"
  EXTRA_PACKAGE_OPTIONS PREFIX Boost
  )

mark_as_advanced(Boost_DIR)

package_on_enabled_script(Boost
  "if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER \"4.8\")
  set(_boost_version \${Boost_MAJOR_VERSION}.\${Boost_MINOR_VERSION})
  if(AKANTU_CORE_CXX11 AND _boost_version VERSION_LESS 1.58 AND _boost_version VERSION_GREATER 1.53)
    add_flags(cxx -DBOOST_SPIRIT_USE_PHOENIX_V3)
  else()
    remove_flags(cxx -DBOOST_SPIRIT_USE_PHOENIX_V3)
  endif()
endif()
")
