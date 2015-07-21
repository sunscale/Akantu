#===============================================================================
# @file   CMakeLists.txt
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date   Wed Jul 9 17:22:12 2014
#
# @brief  CMake file for the python wrapping of akantu
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

#===============================================================================
# Configuration
#===============================================================================
package_get_all_definitions(AKANTU_DEFS)
set(AKA_DEFS "")
foreach (def ${AKANTU_DEFS})
  list(APPEND AKA_DEFS "-D${def}")
endforeach()

set(AKANTU_SWIG_FLAGS -w309,325,401,317,509,503 ${AKA_DEFS})
set(AKANTU_SWIG_OUTDIR ${CMAKE_CURRENT_SOURCE_DIR})
set(AKANTU_SWIG_MODULES swig/akantu.i)

#===============================================================================
# Swig wrapper
#===============================================================================
find_package(SWIG REQUIRED)
find_package(PythonLibs)

package_get_all_include_directories(
  AKANTU_LIBRARY_INCLUDE_DIRS
  )

package_get_all_external_informations(
  AKANTU_EXTERNAL_INCLUDE_DIR
  AKANTU_EXTERNAL_LIBRARIES
)

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/swig
  ${PYTHON_INCLUDE_PATH}
  ${AKANTU_LIBRARY_INCLUDE_DIRS}
  ${PROJECT_BINARY_DIR}/src
  ${AKANTU_EXTERNAL_INCLUDE_DIR}
  )

include(CMakeParseArguments)
function(swig_generate_dependencies _module _depedencies)
  set(_dependencies_script "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/_swig_generate_dependencies.cmake")
  file(WRITE ${_dependencies_script} "
set(_include_directories ${_include_directories})
list(APPEND _include_directories \"./\")

set(_dep)
set(_files_to_process \${_module})
while(_files_to_process)
  list(GET _files_to_process 0 _file)
  list(REMOVE_AT _files_to_process 0)
  file(STRINGS \${_file} _file_content REGEX \"^%include *\\\"(.*)\\\"\")

  set(_includes)
  foreach(_line \${_file_content})
    string(REGEX REPLACE \"^%include *\\\"(.*)\\\"\" \"\\\\1\" _inc \${_line})
    if(_inc)
      list(APPEND _includes \${_inc})
    endif()
  endforeach()

  foreach(_include \${_includes})
    unset(_found)
    foreach(_inc_dir \${_include_directories})
      if(EXISTS \${_inc_dir}/\${_include})
        set(_found \${_inc_dir}/\${_include})
        break()
      endif()
    endforeach()

    if(_found)
      list(APPEND _files_to_process \${_found})
      list(APPEND _dep \${_found})
    endif()
  endforeach()
endwhile()

get_filename_component(_module_we \"\${_module}\" NAME_WE)
set(_dependencies_file \${CMAKE_CURRENT_BINARY_DIR}\${CMAKE_FILES_DIRECTORY}/_swig_\${_module_we}_depends.cmake)
file(WRITE \"\${_dependencies_file}\"
  \"set(_swig_\${_module_we}_depends\")
foreach(_d \${_dep})
  file(APPEND \"\${_dependencies_file}\" \"
  \${_d}\")
endforeach()
file(APPEND \"\${_dependencies_file}\" \"
  )\")
")

  get_directory_property(_include_directories INCLUDE_DIRECTORIES)
  get_filename_component(_module_absolute "${_module}" ABSOLUTE)
  get_filename_component(_module_we "${_module}" NAME_WE)

  set(_dependencies_file ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/_swig_${_module_we}_depends.cmake)

  if(EXISTS ${_dependencies_file})
    include(${_dependencies_file})
  else()
    execute_process(COMMAND ${CMAKE_COMMAND}
      -D_module=${_module_absolute}
      -P ${_dependencies_script}
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    include(${_dependencies_file})
  endif()

  add_custom_command(OUTPUT ${_dependencies_file}
    COMMAND ${CMAKE_COMMAND}
    -D_module=${_module_absolute}
    -P ${_dependencies_script}
    COMMENT "Scannong dependencies for swig module ${_module_we}"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    MAIN_DEPENDENCY ${_module_absolute}
    DEPENDS ${_swig_${_module_we}_depends}
    )


  set(${_depedencies} ${_dependencies_file} PARENT_SCOPE)
endfunction()

function(swig_generate_wrappers project _wrappers)
  cmake_parse_arguments(_swig_opt "" "OUTPUT_DIR" "EXTRA_FLAGS" ${ARGN})

  if(_swig_opt_OUTPUT_DIR)
    set(_output_dir ${_swig_opt_OUTPUT_DIR})
  else()
    set(_output_dir ${CMAKE_CURRENT_BINARY_DIR})
  endif()

  set(_swig_wrappers)
  get_directory_property(_include_directories INCLUDE_DIRECTORIES)
  if(_include_directories)
    string(REPLACE ";" ";-I" _swig_include_directories "${_include_directories}")
  endif()

  foreach(_module ${_swig_opt_UNPARSED_ARGUMENTS})
    swig_generate_dependencies(${_module} _module_dependencies)

    get_filename_component(_module_absolute "${_module}" ABSOLUTE)
    get_filename_component(_module_path "${_module_absolute}" PATH)
    get_filename_component(_module_name "${_module}" NAME)
    get_filename_component(_module_we "${_module}" NAME_WE)
    set(_wrapper "${_output_dir}/${_module_we}_wrapper.cc")
    set(_extra_wrapper "${_output_dir}/${_module_we}.py")
    set(_extra_wrapper_bin "${CMAKE_CURRENT_BINARY_DIR}/${_module_we}.py")


    if(SWIG_FOUND)
      set_source_files_properties("${_wrapper}" PROPERTIES GENERATED 1)
      set_source_files_properties("${_extra_wrapper}" PROPERTIES GENERATED 1)

      set(_dependencies_file ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/_swig_${_module_we}_depends.cmake)

      set(_ouput "${_wrapper}" "${_extra_wrapper}")
      add_custom_command(
	OUTPUT ${_ouput}
	COMMAND "${SWIG_EXECUTABLE}"
	ARGS -python -c++
	${_swig_opt_EXTRA_FLAGS}
	-outdir ${_output_dir}
	-I${_swig_include_directories} -I${_module_path}
	-o "${_wrapper}"
	"${_module_absolute}"
	COMMAND ${CMAKE_COMMAND} -E copy_if_different ${_extra_wrapper} ${_extra_wrapper_bin}
#	MAIN_DEPENDENCY "${_module_absolute}"
        DEPENDS ${_module_dependencies}
	COMMENT "Generating swig wrapper ${_module} -> ${_wrapper}"
	)

      list(APPEND _swig_wrappers ${_wrapper})
    else()
      if(NOT EXISTS ${_wrapper} OR NOT EXISTS "${_extra_wrapper}")
	message(FATAL_ERROR "The file ${_wrapper} and/or ${_extra_wrapper} does not exists and they cannot be generated. Install swig in order to generate them")
      else()
	list(APPEND _swig_wrappers ${_wrapper})
      endif()
    endif()
  endforeach()

  add_custom_target(${project}_generate_swig_wrappers DEPENDS ${_swig_wrappers})

  set(${_wrappers} ${_swig_wrappers} PARENT_SCOPE)
endfunction()


swig_generate_wrappers(akantu AKANTU_SWIG_WRAPPERS
  ${AKANTU_SWIG_MODULES}
  EXTRA_FLAGS ${AKANTU_SWIG_FLAGS})

if(AKANTU_SWIG_WRAPPERS)
  add_library(_akantu MODULE ${AKANTU_SWIG_WRAPPERS})
  target_link_libraries(_akantu akantu ${PYTHON_LIBRARIES})

  set_target_properties(_akantu PROPERTIES PREFIX "")
endif()

