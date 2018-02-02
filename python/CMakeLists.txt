#===============================================================================
# @file   CMakeLists.txt
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Fri Dec 12 2014
# @date last modification: Mon Jan 18 2016
#
# @brief  CMake file for the python wrapping of akantu
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

#===============================================================================
# Configuration
#===============================================================================
package_get_all_definitions(AKANTU_DEFS)
list(REMOVE_ITEM AKANTU_DEFS AKANTU_CORE_CXX11)
#message(${AKANTU_DEFS})
set(AKA_DEFS "")
foreach (def ${AKANTU_DEFS})
  list(APPEND AKA_DEFS "-D${def}")
endforeach()

set(AKANTU_SWIG_FLAGS -w309,325,401,317,509,503,383,384 ${AKA_DEFS})
set(AKANTU_SWIG_OUTDIR ${CMAKE_CURRENT_SOURCE_DIR})
set(AKANTU_SWIG_MODULES swig/akantu.i)

#===============================================================================
# Swig wrapper
#===============================================================================

set(SWIG_REQURIED_VERISON 3.0)
find_package(SWIG ${SWIG_REQURIED_VERISON})
find_package(PythonInterp 2.7 REQUIRED)
mark_as_advanced(SWIG_EXECUTABLE)

package_get_all_include_directories(
  AKANTU_LIBRARY_INCLUDE_DIRS
  )

set(_swig_include_dirs
  ${CMAKE_CURRENT_SOURCE_DIR}/swig
  ${AKANTU_LIBRARY_INCLUDE_DIRS}
  ${PROJECT_BINARY_DIR}/src
  ${AKANTU_EXTERNAL_INCLUDE_DIR}
  )

include(CMakeParseArguments)

function(swig_generate_dependencies _module _depedencies _depedencies_out)
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

  set(${_depedencies_out} ${_swig_${_module_we}_depends} PARENT_SCOPE)

  add_custom_command(OUTPUT ${_dependencies_file}
    COMMAND ${CMAKE_COMMAND}
    -D_module=${_module_absolute}
    -P ${_dependencies_script}
    COMMENT "Scanning dependencies for swig module ${_module_we}"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    MAIN_DEPENDENCY ${_module_absolute}
    DEPENDS ${_swig_${_module_we}_depends}
    )

  set(${_depedencies} ${_dependencies_file} PARENT_SCOPE)
endfunction()

function(swig_generate_wrappers project _wrappers_cpp _wrappers_py)
  cmake_parse_arguments(_swig_opt "" "OUTPUT_DIR;DEPENDENCIES" "EXTRA_FLAGS;INCLUDE_DIRECTORIES" ${ARGN})

  if(_swig_opt_OUTPUT_DIR)
    set(_output_dir ${_swig_opt_OUTPUT_DIR})
  else()
    set(_output_dir ${CMAKE_CURRENT_BINARY_DIR})
  endif()

  set(_swig_wrappers)
  get_directory_property(_include_directories INCLUDE_DIRECTORIES)
  list(APPEND _include_directories ${_swig_opt_INCLUDE_DIRECTORIES})
  if(_include_directories)
    string(REPLACE ";" ";-I" _swig_include_directories "${_include_directories}")
  endif()

  foreach(_module ${_swig_opt_UNPARSED_ARGUMENTS})
    swig_generate_dependencies(${_module} _module_dependencies _depends_out)
    if(_swig_opt_DEPENDENCIES)
      set(${_swig_opt_DEPENDENCIES} ${_depends_out} PARENT_SCOPE)
    endif()

    get_filename_component(_module_absolute "${_module}" ABSOLUTE)
    get_filename_component(_module_path "${_module_absolute}" PATH)
    get_filename_component(_module_name "${_module}" NAME)
    get_filename_component(_module_we "${_module}" NAME_WE)
    set(_wrapper "${_output_dir}/${_module_we}_wrapper.cpp")
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
      list(APPEND _swig_wrappers_py "${_extra_wrapper_bin}")
    else()
      if(NOT EXISTS ${_wrapper} OR NOT EXISTS "${_extra_wrapper}")
        message(FATAL_ERROR "The file ${_wrapper} and/or ${_extra_wrapper} does "
          "not exists and they cannot be generated. Install swig ${SWIG_REQURIED_VERISON} "
          " in order to generate them. Or get them from a different machine, "
          "in order to be able to compile the python interface")
      else()
        list(APPEND _swig_wrappers "${_wrapper}")
        list(APPEND _swig_wrappers_py "${_extra_wrapper_bin}")
      endif()
    endif()
  endforeach()

  add_custom_target(${project}_generate_swig_wrappers DEPENDS ${_swig_wrappers})

  set(${_wrappers_cpp} ${_swig_wrappers} PARENT_SCOPE)
  set(${_wrappers_py} ${_swig_wrappers_py} PARENT_SCOPE)
endfunction()


swig_generate_wrappers(akantu AKANTU_SWIG_WRAPPERS_CPP AKANTU_WRAPPERS_PYTHON
  ${AKANTU_SWIG_MODULES}
  EXTRA_FLAGS ${AKANTU_SWIG_FLAGS}
  DEPENDENCIES _deps
  INCLUDE_DIRECTORIES ${_swig_include_dirs})

if(AKANTU_SWIG_WRAPPERS_CPP)
  string(REPLACE ";" "', '" _ext_files "${AKANTU_SWIG_WRAPPERS_CPP}")
  string(REPLACE ";" "', '" _inc_dirs "${_swig_include_dirs}")
  if (AKANTU_EXTRA_CXX_FLAGS)
    string(STRIP "${AKANTU_EXTRA_CXX_FLAGS}" _tmp_flags)
    string(REGEX REPLACE " +" "', '" _flags "${_tmp_flags} -Wno-maybe-uninitialized")
  else()
    set(_flags "-Wno-maybe-uninitialized")
  endif()

  if(CMAKE_VERBOSE_MAKEFILE)
    set(_quiet)
  else()
    set(_quiet --quiet)
  endif()

  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/setup.py "
from distutils.core import setup
from distutils.core import setup, Extension
import os


os.environ['CC'] = '${CMAKE_CXX_COMPILER}'
os.environ['CXX'] = '${CMAKE_CXX_COMPILER}'
setup(name='akantu',
      license='LGPLv3',
      version='${AKANTU_VERSION}',
      py_modules=['akantu'],
      ext_modules=[Extension('_akantu', ['${_ext_files}'],
                             include_dirs=['${_inc_dirs}'],
                             language='c++',
                             libraries=['akantu${CMAKE_DEBUG_POSTFIX}'],
                             library_dirs=['${PROJECT_BINARY_DIR}/src'],
                             extra_compile_args=['${_flags}', '-std=c++14'])],
      )
")

  add_custom_target(_akantu ALL
    COMMAND ${PYTHON_EXECUTABLE} ./setup.py ${_quiet} --no-user-cfg build_ext --inplace
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    DEPENDS ${AKANTU_SWIG_WRAPPERS_CPP} akantu
    COMMENT "Building akantu's python interface"
    )

  set_directory_properties(PROPERTIES
    ADDITIONAL_MAKE_CLEAN_FILES
      ${CMAKE_CURRENT_BINARY_DIR}/_akantu${CMAKE_SHARED_MODULE_SUFFIX})

  # add_library(_akantu MODULE ${AKANTU_SWIG_WRAPPERS_CPP})
  # target_link_libraries(_akantu akantu)

  # set_target_properties(_akantu PROPERTIES PREFIX "")

  # list(APPEND AKANTU_EXPORT_LIST _akantu)

  # install(TARGETS _akantu
  #   EXPORT ${AKANTU_TARGETS_EXPORT}
  #   LIBRARY DESTINATION lib COMPONENT python NAMELINK_SKIP # for real systems
  #   ARCHIVE DESTINATION lib COMPONENT python
  #   RUNTIME DESTINATION bin COMPONENT python # for windows ...
  #   )

  # install(FILES ${AKANTU_WRAPPERS_PYTHON}
  #   DESTINATION lib COMPONENT python
  #   )

  install(CODE "execute_process(
    COMMAND ${PYTHON_EXECUTABLE} ./setup.py ${_quiet} install --prefix=${AKANTU_PYTHON_INSTALL_PREFIX}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )"
    COMPONENT python)

  package_declare_extra_files_to_package(python_interface ${_deps} ${PROJECT_SOURCE_DIR}/python/${AKANTU_SWIG_MODULES})
endif()
