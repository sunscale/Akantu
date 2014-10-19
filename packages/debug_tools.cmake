#===============================================================================
# @file   debug_tools.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @brief  some usefull tools to debug
#
# @section LICENSE
#
# Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
option(AKANTU_DEBUG_TOOLS "Compiles akantu with the advanced debug tools activated (/!\ can slow down the code)" OFF)
mark_as_advanced(AKANTU_DEBUG_TOOLS)

set(AKANTU_DEBUG_TOOLS_FILES
  common/aka_debug_tools.hh
  )