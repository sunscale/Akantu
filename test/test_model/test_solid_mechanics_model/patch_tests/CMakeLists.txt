#===============================================================================
# @file   CMakeLists.txt
#
# @author David Simon Kammer <david.kammer@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Fri Oct 22 2010
# @date last modification: Mon Dec 07 2015
#
# @brief  configuration for patch tests
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
# @section DESCRIPTION
#
#===============================================================================
set(_meshes)
package_get_variable(ET_ELEMENT_TYPES core _element_list)
foreach(_et ${_element_list})
  if(NOT _et STREQUAL _point_1)
    list(APPEND _meshes data/${_et}.msh)
  endif()
endforeach()

register_gtest_sources(SOURCES patch_test_linear_elastic_explicit.cc
  PACKAGE core
  FILES_TO_COPY data/material_check_stress_plane_strain.dat data/material_check_stress_plane_stress.dat)
register_gtest_sources(SOURCES patch_test_linear_elastic_implicit.cc
  PACKAGE implicit
  FILES_TO_COPY data/material_check_stress_plane_strain.dat data/material_check_stress_plane_stress.dat)
register_gtest_sources(SOURCES patch_test_linear_anisotropic_explicit.cc
  PACKAGE core
  FILES_TO_COPY data/material_anisotropic.dat
  UNSTABLE)
register_gtest_sources(SOURCES test_lumped_mass.cc
  PACKAGE core
  FILES_TO_COPY data/material_lumped.dat)

register_gtest_test(patch_test_linear_elastic
  FILES_TO_COPY ${_meshes})
