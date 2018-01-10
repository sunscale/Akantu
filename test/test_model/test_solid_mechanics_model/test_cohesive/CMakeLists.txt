#===============================================================================
# @file   CMakeLists.txt
#
# @author Marco Vocialta <marco.vocialta@epfl.ch>
#
# @date creation: Fri Oct 22 2010
# @date last modification: Thu Jan 14 2016
#
# @brief  configuration for cohesive elements tests
#
# @section LICENSE
#
# Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
# Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
# Solides)
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

#add_akantu_test(test_cohesive_intrinsic "test_cohesive_intrinsic")
#add_akantu_test(test_cohesive_extrinsic "test_cohesive_extrinsic")
#add_akantu_test(test_cohesive_intrinsic_impl "test_cohesive_intrinsic_impl")
#add_akantu_test(test_cohesive_1d_element "test_cohesive_1d_element")
#add_akantu_test(test_cohesive_extrinsic_implicit "test_cohesive_extrinsic_implicit")

add_akantu_test(test_cohesive_buildfragments "test_cohesive_buildfragments")
add_akantu_test(test_cohesive_insertion "test_cohesive_insertion")
add_akantu_test(test_cohesive_linear_friction "test_cohesive_linear_friction")
#add_akantu_test(test_parallel_cohesive "parallel_cohesive_test")

add_mesh(cohesive_2d_4 data/cohesive_strait_2D.geo 2 1 OUTPUT _cohesive_2d_4.msh)
add_mesh(cohesive_2d_6 data/cohesive_strait_2D.geo 2 2 OUTPUT _cohesive_2d_6.msh)

register_gtest_sources(
  SOURCES test_cohesive.cc
  PACKAGE cohesive_element
  DEPENDS cohesive_2d_4 cohesive_2d_6
  )

register_gtest_test(test_solid_mechanics_model_cohesive)
#===============================================================================
