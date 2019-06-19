#===============================================================================
# @file   CMakeLists.txt
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
#
# @date creation: Fri Sep 03 2010
# @date last modification: Tue Jan 30 2018
#
# @brief  configuration for model tests
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

add_subdirectory(patch_tests)
add_akantu_test(test_common "Tests the common part of the models")
add_akantu_test(test_solid_mechanics_model "Test for the solid mechanics model")
add_akantu_test(test_structural_mechanics_model "Test for the structural mechanics model")
add_akantu_test(test_heat_transfer_model "Test heat transfer model")
add_akantu_test(test_phase_field_model "Test phase field model")
