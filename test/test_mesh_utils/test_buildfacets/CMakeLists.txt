#===============================================================================
# @file   CMakeLists.txt
#
# @author Mauro Corrado <mauro.corrado@epfl.ch>
# @author Marco Vocialta <marco.vocialta@epfl.ch>
#
# @date creation: Tue Dec 07 2010
# @date last modification: Wed Feb 03 2016
#
# @brief  configuration for build facets test
#
# @section LICENSE
#
# Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# Akantu is free  software: you can redistribute it and/or  modify it under the terms  of the  GNU Lesser  General Public  License as published by  the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more details.
#
# You should  have received  a copy  of the GNU  Lesser General  Public License along with Akantu. If not, see <http://www.gnu.org/licenses/>.
#
# @section DESCRIPTION
#
#===============================================================================

register_test(test_buildfacets_triangle_3
  SOURCES test_buildfacets_triangle_3.cc
  FILES_TO_COPY triangle_3.msh
  PACKAGE cohesive_element
  )

register_test(test_buildfacets_triangle_6
  SOURCES test_buildfacets_triangle_6.cc
  FILES_TO_COPY triangle_6.msh
  PACKAGE cohesive_element
  )

register_test(test_buildfacets_quadrangle_4
  SOURCES test_buildfacets_quadrangle_4.cc
  FILES_TO_COPY quadrangle_4.msh
  PACKAGE cohesive_element
  )

register_test(test_buildfacets_quadrangle_8
  SOURCES test_buildfacets_quadrangle_8.cc
  FILES_TO_COPY quadrangle_8.msh
  PACKAGE cohesive_element
  )

register_test(test_buildfacets_mixed2d_linear
  SOURCES test_buildfacets_mixed2d_linear.cc
  FILES_TO_COPY mixed2d_linear.msh
  PACKAGE cohesive_element
  )

register_test(test_buildfacets_mixed2d_quadratic
  SOURCES test_buildfacets_mixed2d_quadratic.cc
  FILES_TO_COPY mixed2d_quadratic.msh
  PACKAGE cohesive_element
  )

register_test(test_buildfacets_tetrahedron_10
  SOURCES test_buildfacets_tetrahedron_10.cc
  FILES_TO_COPY tetrahedron_10.msh
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

register_test(test_buildfacets_pentahedron_6
  SOURCES test_buildfacets_pentahedron_6.cc
  FILES_TO_COPY pentahedron_6.msh
  PACKAGE cohesive_element
  )

register_test(test_buildfacets_pentahedron_15
  SOURCES test_buildfacets_pentahedron_15.cc
  FILES_TO_COPY pentahedron_15.msh
  PACKAGE cohesive_element
  )

register_test(test_buildfacets_mixed3d_linear
  SOURCES test_buildfacets_mixed3d_linear.cc
  FILES_TO_COPY mixed3d_linear.msh
  PACKAGE cohesive_element
  )

register_test(test_buildfacets_mixed3d_quadratic
  SOURCES test_buildfacets_mixed3d_quadratic.cc
  FILES_TO_COPY mixed3d_quadratic.msh
  PACKAGE cohesive_element
  )
