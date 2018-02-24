#===============================================================================
# @file   CMakeLists.txt
#
# @author Mauro Corrado <mauro.corrado@epfl.ch>
#
# @date creation: Fri Sep 03 2010
# @date last modification: Fri Jul 07 2017
#
# @brief  testing the correct behavior of the friction law included in
# the cohesive linear law, in implicit
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

#register_test(test_cohesive_linear_friction
  #SOURCES test_cohesive_friction.cc
  #FILES_TO_COPY mesh_cohesive_friction.msh material.dat
  #PACKAGE cohesive_element
  #)