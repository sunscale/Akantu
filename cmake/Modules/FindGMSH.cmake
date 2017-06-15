#===============================================================================
# @file   FindGMSH.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Mon Dec 08 2014
# @date last modification: Tue Jan 19 2016
#
# @brief  Find gmsh and delacre the add_mesh macro
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

find_program(GMSH gmsh
  DOC "The mesh generetor gmsh")

mark_as_advanced(GMSH)

if (GMSH)
  execute_process(COMMAND ${GMSH} --version
    ERROR_VARIABLE GMSH_VERSION
    ERROR_STRIP_TRAILING_WHITESPACE)
endif()

find_package(PackageHandleStandardArgs)
find_package_handle_standard_args(GMSH DEFAULT_MSG GMSH)
  find_package_handle_standard_args(GMSH
    REQUIRED_VARS GMSH
    VERSION_VAR GMSH_VERSION)

#===============================================================================
function(ADD_MESH MESH_TARGET GEO_FILE DIM ORDER)
  if(GMSH_FOUND)
    set(arguments
      ${MESH_TARGET} ${GEO_FILE} ${DIM} ${ORDER}
      ${ARGN}
      )

    cmake_parse_arguments(ADD_MESH
      ""
      "OUTPUT"
      ""
      ${arguments}
      )

    set(_geo_file ${CMAKE_CURRENT_SOURCE_DIR}/${GEO_FILE})

    set(_r_geo_file "${GEO_FILE}")

    if(ADD_MESH_OUTPUT)
      set(_msh_file ${CMAKE_CURRENT_BINARY_DIR}/${ADD_MESH_OUTPUT})
      set(_r_msh_file "${ADD_MESH_OUTPUT}")
    else(ADD_MESH_OUTPUT)
      get_filename_component(_msh_file "${GEO_FILE}" NAME_WE)
      set(_msh_file ${CMAKE_CURRENT_BINARY_DIR}/${_msh_file}.msh)
      set(_r_msh_file "${_msh_file.msh}")
    endif(ADD_MESH_OUTPUT)

    if(GMSH_VERSION VERSION_LESS 3.0.0)
      set(OPTIMIZE -optimize)
    endif()

    if(EXISTS ${_geo_file})
      add_custom_command(
        OUTPUT ${_msh_file}
        DEPENDS ${_geo_file}
        COMMAND ${GMSH}
        ARGS -${DIM} -order ${ORDER} ${OPTIMIZE} -o ${_msh_file} ${_geo_file} 2>&1 > /dev/null
        COMMENT "Generating the ${DIM}D mesh ${_r_msh_file} (order ${ORDER}) form the geometry ${_r_geo_file}"
        )

      add_custom_target(${MESH_TARGET}
        DEPENDS ${_msh_file})
      set_target_properties(${MESH_TARGET} PROPERTIES RESSOURCES ${_geo_file})
    else(EXISTS ${_geo_file})
      message(WARNING
        "File ${_geo_file} not found for target ${MESH_TARGET}")
    endif(EXISTS ${_geo_file})
  endif(GMSH_FOUND)
endfunction(ADD_MESH)
