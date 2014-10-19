#===============================================================================
# @file   90_iohelper.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Tue Nov 29 2011
# @date last modification: Tue Sep 02 2014
#
# @brief  package description for iohelper
#
# @section LICENSE
#
# Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

option(AKANTU_USE_THIRD_PARTY_IOHELPER "Automatic download of the IOHelper library" ON)
option(AKANTU_USE_IOHELPER "Add IOHelper support in akantu" ON)

mark_as_advanced(AKANTU_USE_IOHELPER)
mark_as_advanced(AKANTU_USE_THIRD_PARTY_IOHELPER)

if (AKANTU_USE_THIRD_PARTY_IOHELPER AND AKANTU_USE_IOHELPER)
  set(IOHELPER_VERSION "1.1")
  set(IOHELPER_GIT     "https://git.epfl.ch/repo/iohelper.git")

  include(ExternalProject)

  ExternalProject_Add(IOHelper
    PREFIX ${PROJECT_BINARY_DIR}/third-party
    GIT_REPOSITORY ${IOHELPER_GIT}
    CMAKE_ARGS <SOURCE_DIR>/
    CMAKE_CACHE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR> -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE} -DCMAKE_CXX_COMPILER:PATH=${CMAKE_CXX_COMPILER}
    BUILD_COMMAND make
    INSTALL_COMMAND make install
    )

    set_third_party_shared_libirary_name(IOHELPER_LIBRARIES iohelper)
    list(APPEND AKANTU_EXTERNAL_LIBRARIES ${IOHELPER_LIBRARIES})
    list(APPEND AKANTU_EXTERNAL_LIB_INCLUDE_DIR ${PROJECT_BINARY_DIR}/third-party/include/iohelper)

    list(APPEND AKANTU_EXTRA_TARGET_DEPENDENCIES IOHelper)
    set(AKANTU_IOHELPER_INCLUDE_DIR ${PROJECT_BINARY_DIR}/third-party/include/iohelper)
    set(AKANTU_IOHELPER ON)
else()
  add_optional_external_package(IOHelper "Add IOHelper support in akantu" ON)
endif()


set(AKANTU_IOHELPER_FILES
  io/dumper/dumper_iohelper.hh
  io/dumper/dumper_iohelper.cc
  io/dumper/dumper_paraview.hh
  io/dumper/dumper_text.cc
  io/dumper/dumper_text.hh
  io/dumper/dumper_paraview.cc
  io/dumper/dumper_homogenizing_field.hh
  io/dumper/dumper_type_traits.hh
  io/dumper/dumper_compute.hh
  io/dumper/dumper_nodal_field.hh
  io/dumper/dumper_quadrature_points_field.hh
  io/dumper/dumper_variable.hh
  io/dumper/dumper_iterator_helper.hh
  io/dumper/dumper_connectivity_field.hh
  io/dumper/dumper_padding_helper.hh
  io/dumper/dumper_elemental_field.hh
  io/dumper/dumper_element_iterator.hh
  io/dumper/dumper_material_internal_field.hh
  io/dumper/dumper_generic_elemental_field.hh
  io/dumper/dumper_generic_elemental_field_tmpl.hh
  )

set(AKANTU_IOHELPER_DOCUMENTATION "
This package activates the IOHelper facilities withing Akantu. This is mandatory if you want to be able to output Paraview files
as well as any Dumper within Akantu.
")
