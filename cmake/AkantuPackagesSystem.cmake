#===============================================================================
# @file   AkantuPackagesSystem.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Sat Jul 18 2015
# @date last modification: Mon Jan 18 2016
#
# @brief  Addition to the PackageSystem specific for Akantu
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
# Element Types
#===============================================================================
#-------------------------------------------------------------------------------
function(package_declare_elements pkg)
  set(_options
    KIND
    ELEMENT_TYPES
    GEOMETRICAL_TYPES
    INTERPOLATION_TYPES
    GEOMETRICAL_SHAPES
    GAUSS_INTEGRATION_TYPES
    INTERPOLATION_TYPES
    INTERPOLATION_KIND
    FE_ENGINE_LISTS
    )

  cmake_parse_arguments(_opt_pkg
    ""
    ""
    "${_options}"
    ${ARGN})

  foreach(_opt ${_options})
    if(_opt_pkg_${_opt})
      package_set_variable(ET_${_opt} ${pkg} ${_opt_pkg_${_opt}})
    endif()
  endforeach()
endfunction()

#-------------------------------------------------------------------------------
function(_transfer_list_to_enum types enum)
  if(${enum})
    set(_enum_tmp ${${enum}})
  else()
    unset(_enum_tmp)
  endif()

  foreach(_type ${${types}})
    # defining the types for the enum
    if(DEFINED _enum_tmp)
      set(_enum_tmp "${_enum_tmp}
  ${_type},")
    else()
      set(_enum_tmp "${_type},")
    endif()
  endforeach()

  set(${enum} ${_enum_tmp} PARENT_SCOPE)
endfunction()


#-------------------------------------------------------------------------------
function(_transfer_list_to_boost_seq types boost_seq)
  if(${boost_seq})
    set(_boost_seq_tmp ${${boost_seq}})
  endif()

  foreach(_type ${${types}})
    if(DEFINED _boost_seq_tmp)
      set(_boost_seq_tmp "${_boost_seq_tmp}${_tabs}\\
  (${_type})")
    else()
      set(_boost_seq_tmp "  (${_type})")
    endif()

    string(LENGTH "${_type}" _length)
    if(_length LESS 13)
      set(_tabs "\t\t\t\t")
    elseif(_length LESS 28)
      set(_tabs "\t\t\t")
    else()
      set(_tabs "\t\t")
    endif()
  endforeach()

  set(${boost_seq} ${_boost_seq_tmp} PARENT_SCOPE)
endfunction()

#-------------------------------------------------------------------------------
function(package_get_element_lists)
  package_get_all_activated_packages(_activated_list)

  set(_lists
    KIND
    ELEMENT_TYPES
    GEOMETRICAL_TYPES
    INTERPOLATION_TYPES
    GEOMETRICAL_SHAPES
    GAUSS_INTEGRATION_TYPES
    INTERPOLATION_TYPES
    INTERPOLATION_KIND
    FE_ENGINE_LISTS
    )

  set(_element_kind "#define AKANTU_ELEMENT_KIND")
  set(_all_element_types "#define AKANTU_ALL_ELEMENT_TYPE\t")

  set(_inter_types_boost_seq "#define AKANTU_INTERPOLATION_TYPES\t\t")

  foreach(_pkg_name ${_activated_list})
    foreach(_list ${_lists})
      string(TOLOWER "${_list}" _l_list)
      _package_get_variable(ET_${_list} ${_pkg_name} _${_l_list})
      _transfer_list_to_enum(_${_l_list} _${_l_list}_enum)
    endforeach()

    if(_interpolation_types)
      _transfer_list_to_boost_seq(_interpolation_types _inter_types_boost_seq)
    endif()

    if(_kind)
      string(TOUPPER "${_kind}" _u_kind)
      if(_element_types)
        set(_boosted_element_types "${_boosted_element_types}
#define AKANTU_ek_${_kind}_ELEMENT_TYPE\t")
        _transfer_list_to_boost_seq(_element_types _boosted_element_types)
        set(_boosted_element_types "${_boosted_element_types}\n")

        # defininf the kinds variables
        set(_element_kinds "${_element_kinds}
#define AKANTU_${_u_kind}_KIND\t(_ek_${_kind})")

        # defining the full list of element
        set(_all_element_types "${_all_element_types}\t\\
  AKANTU_ek_${_kind}_ELEMENT_TYPE")
      endif()

      # defining the full list of kinds
      set(_element_kind "${_element_kind}${_kind_tabs}\t\t\\
  AKANTU_${_u_kind}_KIND")
      set(_kind_tabs "\t")

      # defining the macros
      set(_boost_macros "${_boost_macros}
#define AKANTU_BOOST_${_u_kind}_ELEMENT_SWITCH(macro)                   \\
 AKANTU_BOOST_ELEMENT_SWITCH(macro,                                     \\
                              AKANTU_ek_${_kind}_ELEMENT_TYPE)

#define AKANTU_BOOST_${_u_kind}_ELEMENT_LIST(macro)                     \\
  AKANTU_BOOST_APPLY_ON_LIST(macro,                                     \\
                             AKANTU_ek_${_kind}_ELEMENT_TYPE)
")

      list(APPEND _aka_fe_lists ${_fe_engine_lists})
      foreach(_fe_engine_list ${_fe_engine_lists})
        if(NOT DEFINED _fe_engine_list_${_fe_engine_list})
          string(TOUPPER "${_fe_engine_list}" _u_list)
          string(LENGTH "#define AKANTU_FE_ENGINE_LIST_${_u_list}" _length)
          math(EXPR _length "72 - ${_length}")
          set(_space "")
          while(_length GREATER 0)
            if(CMAKE_VERSION VERSION_GREATER 3.0)
              string(CONCAT _space "${_space}" " ")
            else()
              set(_space "${_space} ")
            endif()
            math(EXPR _length "${_length} - 1")
          endwhile()

          set(_fe_engine_list_${_fe_engine_list}
            "#define AKANTU_FE_ENGINE_LIST_${_u_list}${_space}\\
  AKANTU_GENERATE_KIND_LIST((_ek_${_kind})")
        else()
          set(_fe_engine_list_${_fe_engine_list}
            "${_fe_engine_list_${_fe_engine_list}}\t\t\t\t\\
                            (_ek_${_kind})")
        endif()
      endforeach()
    endif()
  endforeach()

  if(_aka_fe_lists)
    list(REMOVE_DUPLICATES _aka_fe_lists)
    foreach(_fe_list ${_aka_fe_lists})
      set(_aka_fe_defines "${_fe_engine_list_${_fe_list}})\n${_aka_fe_defines}")
    endforeach()
  endif()

  foreach(_list ${_lists})
    string(TOLOWER "${_list}" _l_list)
    set(AKANTU_${_list}_ENUM ${_${_l_list}_enum} PARENT_SCOPE)
  endforeach()

  set(AKANTU_INTERPOLATION_TYPES_BOOST_SEQ ${_inter_types_boost_seq} PARENT_SCOPE)
  set(AKANTU_ELEMENT_TYPES_BOOST_SEQ ${_boosted_element_types} PARENT_SCOPE)
  set(AKANTU_ELEMENT_KINDS_BOOST_SEQ ${_element_kinds} PARENT_SCOPE)
  set(AKANTU_ELEMENT_KIND_BOOST_SEQ ${_element_kind} PARENT_SCOPE)
  set(AKANTU_ALL_ELEMENT_BOOST_SEQ ${_all_element_types} PARENT_SCOPE)
  set(AKANTU_ELEMENT_KINDS_BOOST_MACROS ${_boost_macros} PARENT_SCOPE)
  set(AKANTU_FE_ENGINE_LISTS ${_aka_fe_defines} PARENT_SCOPE)
endfunction()

#-------------------------------------------------------------------------------
function(package_get_element_types pkg list)
  package_get_name(${pkg} _pkg_name)
  _package_get_variable(ET_ELEMENT_TYPES ${_pkg_name} _tmp_list)
  set(${list} ${_tmp_list} PARENT_SCOPE)
endfunction()

#===============================================================================
# Material specific
#===============================================================================
#-------------------------------------------------------------------------------
function(package_declare_material_infos pkg)
  cmake_parse_arguments(_opt_pkg
    ""
    ""
    "LIST;INCLUDE"
    ${ARGN})

  package_set_variable(MATERIAL_LIST ${pkg} ${_opt_pkg_LIST})
  package_set_variable(MATERIAL_INCLUDE ${pkg} ${_opt_pkg_INCLUDE})
endfunction()

#-------------------------------------------------------------------------------
function(package_get_all_material_includes includes)
  _package_get_variable_for_activated(MATERIAL_INCLUDE _includes)

  foreach(_mat_inc ${_includes})
    if(DEFINED _mat_includes)
      set(_mat_includes "${_mat_includes}\n#include \"${_mat_inc}\"")
    else()
      set(_mat_includes "#include \"${_mat_inc}\"")
    endif()
  endforeach()

  set(${includes} ${_mat_includes} PARENT_SCOPE)
endfunction()

#-------------------------------------------------------------------------------
function(package_get_all_material_lists lists)
  _package_get_variable_for_activated(MATERIAL_LIST _lists)

  foreach(_mat_list ${_lists})
    if(DEFINED _mat_lists)
      set(_mat_lists "${_mat_lists}\n  ${_mat_list}\t\t\t\\")
    else()
      set(_mat_lists "  ${_mat_list}\t\t\t\\")
    endif()
  endforeach()

  set(${lists} ${_mat_lists} PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Extra files to consider in source package generated by CPack
# ------------------------------------------------------------------------------
function(package_declare_extra_files_to_package pkg)
  set(_types SOURCES MANUAL TESTS PROJECT)
  cmake_parse_arguments(_extra_files
    ""
    ""
    "${_types}"
    ${ARGN})

  set(_files ${_extra_files_UNPARSED_ARGUMENTS})

  package_get_sources_folder(${pkg} _folder_SOURCES)
  package_get_manual_folder(${pkg} _folder_MANUAL)
  package_get_tests_folder(${pkg} _folder_TESTS)
  set(_folder_PROJECT ${PROJECT_SOURCE_DIR})

  foreach(_type ${_types})
    if(_extra_files_${_type})
      foreach(_file ${_extra_files_${_type}})
        list(APPEND _files ${_folder_${_type}}/${_file})
        if(NOT EXISTS ${_folder_${_type}}/${_file})
          message(SEND_ERROR "The package ${pkg} tries to register the file ${_file} (as a ${_type} file).
This file cannot be found.")
        endif()
      endforeach()
    endif()
  endforeach()

  package_set_variable(EXTRA_FILES ${pkg} ${_files})
endfunction()

# ------------------------------------------------------------------------------
function(package_add_files_to_package)
  set(_files)
  foreach(_file ${ARGN})
    list(APPEND _files ${PROJECT_SOURCE_DIR}/${_file})
  endforeach()
  package_add_to_project_variable(EXTRA_FILES ${_files})
endfunction()

function(package_get_files_for_package files)
  package_get_project_variable(EXTRA_FILES _tmp)
  set(${files} ${_tmp} PARENT_SCOPE)
endfunction()


package_add_files_to_package(
  .clang-format
  AUTHORS
  README
  VERSION
  COPYING
  COPYING.lesser
  CTestConfig.cmake
  cmake/akantu_environement.sh.in
  cmake/akantu_environement.csh.in
  cmake/akantu_install_environement.sh.in
  cmake/akantu_install_environement.csh.in
  cmake/Modules/CMakeFlagsHandling.cmake
  cmake/Modules/CMakePackagesSystem.cmake
  cmake/Modules/CMakePackagesSystemGlobalFunctions.cmake
  cmake/Modules/CMakePackagesSystemPrivateFunctions.cmake
  cmake/Modules/CMakeVersionGenerator.cmake
  cmake/Modules/PCHgcc.cmake
  cmake/AkantuBuildTreeSettings.cmake.in
  cmake/AkantuConfig.cmake.in
  cmake/AkantuConfigVersion.cmake.in
  cmake/AkantuCPack.cmake
  cmake/AkantuCPackMacros.cmake
  cmake/AkantuInstall.cmake
  cmake/AkantuMacros.cmake
  cmake/AkantuPackagesSystem.cmake
  cmake/AkantuUse.cmake
  cmake/AkantuSimulationMacros.cmake
  cmake/material_lister.cc
  cmake/Modules/FindGMSH.cmake
  )
