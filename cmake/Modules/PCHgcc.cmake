#===============================================================================
# @file   PCHgcc.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Wed Aug 31 2011
# @date last modification: Sun Jan 06 2013
#
# @brief  
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

# (ADD_PCH_RULE  _header_filename _src_list)
# Version 7/26/2010 4:55pm
#
# use this macro before "add_executable"
#
# _header_filename
#	header to make a .gch
#
# _src_list
#   the variable name (do not use ${..}) which contains a
#     a list of sources (a.cpp b.cpp c.cpp ...)
#  This macro will append a header file to it, then this src_list can be used in
#	"add_executable..."
#
#
# Now a .gch file should be generated and gcc should use it.
#		(add -Winvalid-pch to the cpp flags to verify)
#
# make clean should delete the pch file
#
# example : ADD_PCH_RULE(headers.h myprog_SRCS)


function(ADD_PCH_RULE _header_filename _src_list)
  set(_gch_filename "${CMAKE_CURRENT_BINARY_DIR}/${_header_filename}.gch")

  list(APPEND ${_src_list} ${_gch_filename})

  set(_args ${CMAKE_CXX_FLAGS} -D__aka_inline__=inline)

  get_filename_component(_gch_filename_path ${_gch_filename} PATH)
  file(MAKE_DIRECTORY ${_gch_filename_path})

  #  list(APPEND _args -c ${CMAKE_CURRENT_SOURCE_DIR}/${_header_filename})
  list(APPEND _args -c ${CMAKE_CURRENT_SOURCE_DIR}/${_header_filename} -o ${_gch_filename} -Winvalid-pch)

  get_directory_property(DIRINC INCLUDE_DIRECTORIES)
  foreach (_inc ${DIRINC})
    list(APPEND _args "-I" ${_inc})
  endforeach(_inc ${DIRINC})
  separate_arguments(_args)

  set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/${_header_filename} PROPERTIES GENERATED 1)
  set_source_files_properties(${_gch_filename} PROPERTIES GENERATED 1)
  add_custom_command(OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/${_header_filename}
    COMMAND ${CMAKE_COMMAND} -E copy  ${CMAKE_CURRENT_SOURCE_DIR}/${_header_filename} ${CMAKE_CURRENT_BINARY_DIR}/${_header_filename} # ensure same directory! Required by gcc
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${_header_filename}
    )

  add_custom_command(OUTPUT ${_gch_filename}
    COMMAND ${CMAKE_CXX_COMPILER} ${CMAKE_CXX_COMPILER_ARG1} ${_args}
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${_header_filename}
    )
endfunction(ADD_PCH_RULE _header_filename _src_list)
