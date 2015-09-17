#===============================================================================
# @file   CMakeLists.txt
#
# @author Marco Vocialta <marco.vocialta@epfl.ch>
# @author Mauro Corrado <mauro.corrado@epfl.ch>
# @date   Thu Sep 17 13:01:18 2015
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

register_test(test_buildfacets_tetrahedron
  SOURCES test_buildfacets_tetrahedron.cc
  FILES_TO_COPY tetrahedron.msh
  PACKAGE cohesive_element
  )
register_test(test_buildfacets_hexahedron_8
  SOURCES test_buildfacets_hexahedron_8.cc
  FILES_TO_COPY hexahedron_8.msh
  PACKAGE cohesive_element
  )
register_test(test_buildfacets_hexahedron_20
  SOURCES test_buildfacets_hexahedron_20.cc
  FILES_TO_COPY hexahedron_20.msh
  PACKAGE cohesive_element
  )

