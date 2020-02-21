#===============================================================================
# @file   mpi.cmake
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Mon Nov 21 2011
# @date last modification: Wed Jan 20 2016
#
# @brief  package description for mpi
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
set(MPI_CXX_SKIP_MPICXX TRUE)
package_declare(MPI EXTERNAL
  DESCRIPTION "Add MPI support in akantu"
  EXTRA_PACKAGE_OPTIONS PREFIX MPI_C MPI
  )

package_declare_sources(MPI
  synchronizer/mpi_communicator_data.hh
  synchronizer/communicator_mpi_inline_impl.cc
  )

function(_add_to_mpi_preflags flag)
  if(NOT MPIEXEC_PREFLAGS MATCHES "${flag}")
	  string(STRIP "${flag} ${MPIEXEC_PREFLAGS}" _preflags)
	  set(MPIEXEC_PREFLAGS "${_preflags}" CACHE STRING "" FORCE)
	endif()
endfunction()

function(add_extra_mpi_options)
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
        set(MPI_STD_VERSION "${_mpi_VERSION}.${_mpi_SUBVERSION}" CACHE INTERNAL "")
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
          set(MPI_ID_VERSION "${MPI_VERSION_MAJOR}.${MPI_VERSION_MINOR}.${MPI_VERSION_RELEASE}"
	          CACHE INTERNAL "")
          _add_to_mpi_preflags("--oversubscribe")
          if(AKANTU_RUN_IN_DOCKER)
	          _add_to_mpi_preflags("--allow-run-as-root")
	        endif()
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


    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      set( _flags "${_flags} -Wno-literal-suffix")
    endif()
  endif()

  include(FindPackageMessage)
  if(MPI_FOUND)
    find_package_message(MPI "MPI ID: ${MPI_ID} ${MPI_ID_VERSION} (MPI standard ${MPI_STD_VERSION})" "${MPI_STD_VERSION}")
  endif()

  set(MPI_EXTRA_COMPILE_FLAGS "${_flags}" CACHE STRING "Extra flags for MPI" FORCE)
  mark_as_advanced(MPI_EXTRA_COMPILE_FLAGS)

  #package_get_source_files(MPI _srcs _pub _priv)
  #list(APPEND _srcs "common/aka_error.cc")

  #set_property(SOURCE ${_srcs} PROPERTY COMPILE_FLAGS "${_flags}")
  package_set_compile_flags(MPI CXX ${_flags})
endfunction()

package_on_enabled_script(MPI
  "
add_extra_mpi_options()
mask_package_options(MPI)
"
)

package_declare_documentation(MPI
  "This is a meta package providing access to MPI."
  ""
  "Under Ubuntu (14.04 LTS) the installation can be performed using the commands:"
  "\\begin{command}"
  "  > sudo apt-get install libopenmpi-dev"
  "\\end{command}"
  ""
  "Under Mac OS X the installation requires the following steps:"
  "\\begin{command}"
  "  > sudo port install mpich-devel"
  "\\end{command}"
  )

package_set_package_system_dependency(MPI deb mpi-default-bin)
package_set_package_system_dependency(MPI deb-src mpi-default-dev)
