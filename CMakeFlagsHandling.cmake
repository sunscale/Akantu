#===============================================================================
# @file   CMakeFlagsHandling.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date   Tue Oct 16 14:05:02 2012
#
# @brief  Set of macros used by akantu to handle the package system
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

if(_CMAKE_FLAGS_HANDLING)
  return()
endif()
set(_CMAKE_FLAGS_HANDLING TRUE)


#===============================================================================
# Compilation options handling
#===============================================================================
macro(_get_flags_message type desc)
  if(${type} MATCHES "C..")
    set(${desc} "Flags used by the compiler during all build types.")
  elseif(${type} MATCHES "EXE_LINKER")
    set(${desc} "Flags used by the linker.")
  elseif(${type} MATCHES "SHARED_LINKER")
    set(${desc} "Flags used by the linker during the creation of dll's.")
  endif()
endmacro()

#===============================================================================
macro(add_flags type flag)
  string(TOUPPER ${type} _type)
  set(_var CMAKE_${_type}_FLAGS)
  _get_flags_message(${_type} _desc)
  if(NOT ${_var} MATCHES "${flag}")
    set(${_var} "${flag} ${${_var}}" CACHE STRING ${_desc} FORCE)
  endif()
endmacro()

#===============================================================================
macro(remove_flags type flag)
  string(TOUPPER ${type} _type)
  set(_var CMAKE_${_type}_FLAGS)
  _get_flags_message(${_type} _desc)
  string(REPLACE "${flag} " "" ${_var} "${${_var}}")
  set(${_var} "${${_var}}" CACHE STRING ${_desc} FORCE)
endmacro()
#===============================================================================
