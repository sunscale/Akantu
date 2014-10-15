#===============================================================================
# @file   CMakeLists.txt
#
# @author Marco Vocialta <marco.vocialta@epfl.ch>
#
# @date   Tue Oct 28 09:43:18 2013
#
# @brief  configuration for parallel test of facet stress synchronizer
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
add_mesh(test_cohesive_facet_stress_synchronizer_mesh tetrahedron.geo 3 2)

if(AKANTU_PARALLEL)
  register_test(test_cohesive_facet_stress_synchronizer
    SOURCES test_cohesive_facet_stress_synchronizer.cc
    DEPENDENCIES test_cohesive_facet_stress_synchronizer_mesh
    PACKAGE cohesive_element
    FILES_TO_COPY material.dat)
endif()
