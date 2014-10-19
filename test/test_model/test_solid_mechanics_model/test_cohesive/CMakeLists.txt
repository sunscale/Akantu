#===============================================================================
# @file   CMakeLists.txt
#
# @author Marco Vocialta <marco.vocialta@epfl.ch>
#
# @date creation: Tue May 08 2012
# @date last modification: Fri Jun 21 2013
#
# @brief  configuration for cohesive elements tests
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
# @section DESCRIPTION
#
#===============================================================================

add_akantu_test(test_cohesive_buildfacets "test_cohesive_buildfacets" PACKAGE cohesive_element)
add_akantu_test(test_cohesive_intrinsic "test_cohesive_intrinsic" PACKAGE cohesive_element)
add_akantu_test(test_cohesive_extrinsic "test_cohesive_extrinsic" PACKAGE cohesive_element)
add_akantu_test(test_cohesive_buildfragments "test_cohesive_buildfragments" PACKAGE cohesive_element)
add_akantu_test(test_cohesive_intrinsic_impl "test_cohesive_intrinsic_impl" PACKAGE cohesive_element)
add_akantu_test(test_cohesive_1d_element "test_cohesive_1d_element" PACKAGE cohesive_element)
#===============================================================================
