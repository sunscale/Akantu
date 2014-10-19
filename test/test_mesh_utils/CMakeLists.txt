#===============================================================================
# @file   CMakeLists.txt
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Mon Feb 07 2011
# @date last modification: Tue Nov 06 2012
#
# @brief  configuration for MeshUtils tests
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
# @section DESCRIPTION
#
#===============================================================================

#===============================================================================
# List of tests
#===============================================================================
add_akantu_test(test_mesh_io "Test mesh io object")
add_akantu_test(test_facet_extraction "Test mesh utils facet extraction")
add_akantu_test(test_pbc_tweak "Test pbc facilities")

register_test(test_purify_mesh
  SOURCES test_purify_mesh.cc
  FILES_TO_COPY purify_mesh.msh DESTINATION
  )

add_akantu_test(test_mesh_partitionate "Test mesh partition creation" PACKAGE scotch)


