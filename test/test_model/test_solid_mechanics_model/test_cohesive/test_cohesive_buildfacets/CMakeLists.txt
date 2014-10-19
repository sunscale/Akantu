#===============================================================================
# @file   CMakeLists.txt
#
# @author Marco Vocialta <marco.vocialta@epfl.ch>
#
# @date   Tue May 08 13:01:18 2012
#
# @brief  configuration for build facets test
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
# @section DESCRIPTION
#
#===============================================================================

register_test(test_cohesive_buildfacets_tetrahedron
  SOURCES test_cohesive_buildfacets_tetrahedron.cc
  FILES_TO_COPY tetrahedron.msh)
register_test(test_cohesive_buildfacets_hexahedron
  SOURCES test_cohesive_buildfacets_hexahedron.cc
  FILES_TO_COPY hexahedron.msh)

#add_mesh(test_cohesive_buildfacets_mesh mesh.geo 3 2)
#add_dependencies(test_cohesive_buildfacets test_cohesive_buildfacets_mesh)

