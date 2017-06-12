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

add_mesh(test_weight_hexa_mesh colone_hexa.geo 3 1)
add_mesh(test_weight_tetra_mesh colone_tetra.geo 3 1)

register_test(test_weight
  SOURCES colone_weight.cc
  DEPENDS test_weight_hexa_mesh test_weight_tetra_mesh
  FILES_TO_COPY material_colone.dat
  DIRECTORIES_TO_CREATE paraview/test_weight_hexa paraview/test_weight_tetra
  PACKAGE core
  )

#===============================================================================
set(LIST_TYPES
  segment_2
  segment_3
  tetrahedron_4
  tetrahedron_10
  hexahedron_8
  hexahedron_20
  pentahedron_6
  pentahedron_15
  )

set(LIST_TYPES_2D
  triangle_3
  triangle_6
  quadrangle_4
  quadrangle_8
)

add_akantu_test(explicit "Explicit patch tests")
add_akantu_test(implicit "Implicit patch tests")

add_akantu_test(lumped_mass "Test the mass lumping")
