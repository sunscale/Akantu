#===============================================================================
# @file   CMakePackagesSystem.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @brief  Addition to the PackageSystem specific for Akantu
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
#define AKANTU_BOOST_${_u_kind}_ELEMENT_SWITCH(macro)			\\
 AKANTU_BOOST_ELEMENT_SWITCH(macro,					\\
			      AKANTU_ek_${_kind}_ELEMENT_TYPE)

#define AKANTU_BOOST_${_u_kind}_ELEMENT_LIST(macro)			\\
  AKANTU_BOOST_APPLY_ON_LIST(macro,					\\
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

#   package_get_all_deactivated_packages(_deactivated_list)
#   foreach(_pkg_name ${_deactivated_list})
#     _package_get_variable(ET_KIND ${_pkg_name} _kind)
#     if(_kind)
#       string(TOUPPER "${_kind}" _u_kind)
#       set(_element_kinds "${_element_kinds}
# #define AKANTU_${_u_kind}_KIND")
#     endif()
#   endforeach()

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

#-------------------------------------------------------------------------------
function(add_extra_mpi_options)
  package_is_activated(MPI _act)
  if(_act)
    unset(MPI_ID CACHE)
    package_get_include_dir(MPI _include_dir)
    foreach(_inc_dir ${_include_dir})
      if(EXISTS "${_inc_dir}/mpi.h")
        if(NOT MPI_ID)
          file(STRINGS "${_inc_dir}/mpi.h" _mpi_version REGEX "#define MPI_(SUB)?VERSION .*")
          foreach(_ver ${_mpi_version})
            string(REGEX MATCH "MPI_(VERSION|SUBVERSION) *([0-9]+)" _tmp "${_ver}")
            set(_mpi_${CMAKE_MATCH_1} ${CMAKE_MATCH_2})
          endforeach()
          set(MPI_VERSION "${_mpi_VERSION}.${_mpi_SUBVERSION}" CACHE INTERNAL "")
        endif()

        if(NOT MPI_ID)
          # check if openmpi
          file(STRINGS "${_inc_dir}/mpi.h" _ompi_version REGEX "#define OMPI_.*_VERSION .*")
          if(_ompi_version)
            set(MPI_ID "OpenMPI" CACHE INTERNAL "")
            foreach(_version ${_ompi_version})
              string(REGEX MATCH "OMPI_(.*)_VERSION (.*)" _tmp "${_version}")
              if(_tmp)
                set(MPI_VERSION_${CMAKE_MATCH_1} ${CMAKE_MATCH_2})
              endif()
            endforeach()
            set(MPI_ID_VERSION "${MPI_VERSION_MAJOR}.${MPI_VERSION_MINOR}.${MPI_VERSION_RELEASE}" CACHE INTERNAL "")
          endif()
        endif()

        if(NOT MPI_ID)
          # check if intelmpi
          file(STRINGS "${_inc_dir}/mpi.h" _impi_version REGEX "#define I_MPI_VERSION .*")
          if(_impi_version)
            set(MPI_ID "IntelMPI" CACHE INTERNAL "")
            string(REGEX MATCH "I_MPI_VERSION \"(.*)\"" _tmp "${_impi_version}")
            if(_tmp)
              set(MPI_ID_VERSION "${CMAKE_MATCH_1}" CACHE INTERNAL "")
            endif()
          endif()
        endif()

        if(NOT MPI_ID)
          # check if mvapich2
          file(STRINGS "${_inc_dir}/mpi.h" _mvapich2_version REGEX "#define MVAPICH2_VERSION .*")
          if(_mvapich2_version)
            set(MPI_ID "MPVAPICH2" CACHE INTERNAL "")
            string(REGEX MATCH "MVAPICH2_VERSION \"(.*)\"" _tmp "${_mvapich2_version}")
            if(_tmp)
              set(MPI_ID_VERSION "${CMAKE_MATCH_1}" CACHE INTERNAL "")
            endif()
          endif()
        endif()

        if(NOT MPI_ID)
          # check if mpich (mpich as to be checked after all the mpi that derives from it)
          file(STRINGS "${_inc_dir}/mpi.h" _mpich_version REGEX "#define MPICH_VERSION .*")
          if(_mpich_version)
            set(MPI_ID "MPICH" CACHE INTERNAL "")
            string(REGEX MATCH "I_MPI_VERSION \"(.*)\"" _tmp "${_mpich_version}")
            if(_tmp)
              set(MPI_ID_VERSION "${CMAKE_MATCH_1}" CACHE INTERNAL "")
            endif()
          endif()
        endif()
      endif()
    endforeach()

    if(MPI_ID STREQUAL "IntelMPI" OR
        MPI_ID STREQUAL "MPICH" OR
        MPI_ID STREQUAL "MVAPICH2")
      set(_flags "-DMPICH_IGNORE_CXX_SEEK")
    elseif(MPI_ID STREQUAL "OpenMPI")
      set(_flags "-DOMPI_SKIP_MPICXX")

      package_is_activated(core_cxx11 _act)
      if(_act)
        set( _flags "${_flags} -Wno-literal-suffix")
      endif()
    endif()

    message(STATUS "MPI ID: ${MPI_ID} ${MPI_ID_VERSION} - (MPI spec v${MPI_VERSION})")

    set(MPI_EXTRA_COMPILE_FLAGS "${_flags}" CACHE STRING "Extra flags for MPI" FORCE)
    mark_as_advanced(MPI_EXTRA_COMPILE_FLAGS)

    package_get_source_files(MPI _srcs _pub _priv)
    list(APPEND _srcs "common/aka_error.cc")

    set_property(SOURCE ${_srcs} PROPERTY COMPILE_FLAGS "${_flags}")
  endif()

endfunction()
