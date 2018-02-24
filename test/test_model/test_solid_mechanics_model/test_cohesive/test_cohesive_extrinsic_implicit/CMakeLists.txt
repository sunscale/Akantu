#===============================================================================
# @file   CMakeLists.txt
#
# @author Mauro Corrado <mauro.corrado@epfl.ch>
#
# @date creation: Fri Sep 03 2010
# @date last modification: Wed Feb 03 2016
#
# @brief  checking correct assembling of stiffness matrix in case of cohe elems
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

register_test(test_assembling_K_cohe_elements
  SOURCES test_assembling_K_cohe_elements.cc
  FILES_TO_COPY quadrangle.msh material.dat K_matrix_verified.dat
  PACKAGE cohesive_element
  )
