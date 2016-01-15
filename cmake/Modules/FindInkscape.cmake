#===============================================================================
# @file   FindInkscape.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Sun Oct 19 2014
#
# @brief  find_package for inkscape
#
# @section LICENSE
#
# Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory
# (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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

# Module thaat checks for inkscape
#
# Sets the following variables
#
# INSCAPE: Path to inkscape to generate .png's form .svg's
#
# Provides the following functions:
#
# inkscape_generate_png_from_svg([OUTPUT_DIR <output_dir>] <pngfile1.png> [<pngfile2.png> ....])
#
# Generates pngfile1, ... from svg input files pngfile1.svg, ....
# The output directory can be specified with the option OUTPUT_DIR. If it is omitted
# the files will be generated in CMAKE_CURRENT_BINARY_DIR.

find_program(INKSCAPE inkscape DOC "Path to inkscape to generate png files from svg files")
find_program(CONVERT convert DOC "Path to convert program")
if(INKSCAPE)
  set(INKSCAPE_FOUND True)
endif(INKSCAPE)

include(CMakeParseArguments)

function(inkscape_generate_png_from_svg)
  if(NOT INKSCAPE)
    return()
  endif(NOT INKSCAPE)
  cmake_parse_arguments(INKSCAPE "" "DPI" "" ${ARGN})
  if(NOT INKSCAPE_DPI)
    set(INKSCAPE_DPI 90)
  endif(NOT INKSCAPE_DPI)

  foreach(pic ${INKSCAPE_UNPARSED_ARGUMENTS})
    string(REGEX REPLACE "\\.[a-zA-Z]+" ".png" output ${pic})
    add_custom_command(OUTPUT ${output} 
      COMMAND ${INKSCAPE} --export-dpi=${INKSCAPE_DPI} -e ${output} -f ${pic}
      DEPENDS ${pic} 
      COMMENT "Generating ${output} from ${pic}"
      )
  endforeach(pic)  
endfunction(inkscape_generate_png_from_svg)

function(inkscape_generate_eps_from_svg)
  cmake_parse_arguments(INKSCAPE "" "INPUT_DIR;OUTPUT_DIR;DPI" "" ${ARGN})
  if(NOT INKSCAPE_INPUT_DIR)
    set(INKSCAPE_INPUT_DIR ${CMAKE_CURRENT_SOURCE_DIR})
  endif(NOT INKSCAPE_INPUT_DIR)
  if(NOT INKSCAPE_INPUT_DIR)
    set(INKSCAPE_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR})
  endif(NOT INKSCAPE_INPUT_DIR)

  foreach(_pic ${INKSCAPE_UNPARSED_ARGUMENTS})
    string(REGEX REPLACE "\\.[a-zA-Z]+" ".png" input "${_pic}")
    string(REGEX REPLACE "\\.[a-zA-Z]+" ".svg" svginput "${_pic}")
    
    add_custom_target(${input}
      COMMAND ${INKSCAPE} --export-dpi=${INKSCAPE_DPI} -e ${input} ${CMAKE_CURRENT_SOURCE_DIR}/${svginput}
      COMMENT "Generating ${INKSCAPE_OUTPUT_DIR}/${pic} from ${CMAKE_CURRENT_SOURCE_DIR}/${input}")
    add_custom_command(OUTPUT ${_pic} 
      COMMAND ${CONVERT} ${INKSCAPE_INPUT_DIR}/${input} EPS:${_pic}
      DEPENDS ${input} 
      COMMENT "Converting ${INKSCAPE_INPUT_DIR}/${input} to ${INKSCAPE_OUTPUT_DIR}/${_pic}"
      WORKING_DIRECTORY  ${INKSCAPE_OUTPUT_DIR})
  endforeach(_pic ${INKSCAPE_UNPARSED_ARGUMENTS})
endfunction(inkscape_generate_eps_from_svg)