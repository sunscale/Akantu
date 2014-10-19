#===============================================================================
# @file   CMakeLists.txt
#
# @author Fabian Barras <fabian.barras@epfl.ch>
#
# @date creation: Fri Jul 15 2011
# @date last modification: Mon Jul 07 2014
#
# @brief  Structural Mechanics Model Tests
#
# @section LICENSE
#
# Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#===============================================================================

register_test(test_structural_mechanics_model_bernoulli_beam_2_exemple_1_1
  test_structural_mechanics_model_bernoulli_beam_2_exemple_1_1.cc)

register_test(test_structural_mechanics_model_bernoulli_beam_2_exemple_1_1_y
  test_structural_mechanics_model_bernoulli_beam_2_exemple_1_1_y.cc)

register_test(test_structural_mechanics_model_bernoulli_beam_3_exemple_1_1_xy
  test_structural_mechanics_model_bernoulli_beam_3_exemple_1_1_xy.cc)

register_test(test_structural_mechanics_model_bernoulli_beam_3_exemple_1_1_zy
  test_structural_mechanics_model_bernoulli_beam_3_exemple_1_1_zy.cc)

register_test(test_structural_mechanics_model_bernoulli_beam_3_exercice_12_10_13
  test_structural_mechanics_model_bernoulli_beam_3_exercice_12_10_13.cc)

register_test(test_structural_mechanics_model_bernoulli_beam_dynamics
  test_structural_mechanics_model_bernoulli_beam_dynamics.cc)

register_test(test_structural_mechanics_model_kirchhoff_shell_patch_test_4_5_5
  test_structural_mechanics_model_kirchhoff_shell_patch_test_4_5_5.cc)

add_mesh(complicated_mesh complicated.geo 1 1)
register_test(test_structural_mechanics_model_bernoulli_beam_2_complicated
  test_structural_mechanics_model_bernoulli_beam_2_complicated.cc
  DEPENDENCIES complicated_mesh)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview)
