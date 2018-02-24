#===============================================================================
# @file   CMakeLists.txt
#
# @author Marco Vocialta <marco.vocialta@epfl.ch>
#
# @date creation: Fri Sep 03 2010
# @date last modification: Wed Feb 03 2016
#
# @brief  configuration for parallel test for extrinsic cohesive elements
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

if(AKANTU_PARALLEL)
  add_mesh(test_cohesive_1d_element_mesh bar.geo 1 2)

  register_test(test_cohesive_1d_element
    SOURCES test_cohesive_1d_element.cc
    DEPENDS test_cohesive_1d_element_mesh
    FILES_TO_COPY material.dat
    DIRECTORIES_TO_CREATE paraview
    PACKAGE cohesive_element)
endif()
