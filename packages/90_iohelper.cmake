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

package_declare(IOHelper EXTERNAL
  DESCRIPTION "Add IOHelper support in akantu"
  SYSTEM OFF
  DEFAULT ON)

package_get_name(IOHelper _pkg_name)
package_use_system(${_pkg_name} _use_system)

if(NOT ${_use_system})
  package_get_option_name(${_pkg_name} _option_name)

  if(${_option_name})
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
    package_set_libraries(${_pkg_name} ${IOHELPER_LIBRARIES})
    package_set_include_dir(${_pkg_name} ${PROJECT_BINARY_DIR}/third-party/include/iohelper)

    package_add_extra_dependency(IOHelper IOHelper)
  endif()
endif()

package_declare_sources(IOHelper
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
