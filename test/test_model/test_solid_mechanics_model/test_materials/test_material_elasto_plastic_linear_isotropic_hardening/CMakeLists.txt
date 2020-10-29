#===============================================================================
# @file   CMakeLists.txt
#
# @author Jaehyun Cho <jaehyun.cho@epfl.ch>
#
# @date creation: Fri Oct 22 2010
# @date last modification: Wed Mar 09 2016
#
# @brief  test for material type elasto plastic linear isotropic hardening using
# tension-compression test
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

#===============================================================================
add_mesh(test_material_standard_linear_isotropic_hardening_mesh
  test_material_elasto_plastic_linear_isotropic_hardening.geo 2 1)

register_test(test_material_standard_linear_isotropic_hardening
  SOURCES test_material_elasto_plastic_linear_isotropic_hardening.cc
  DEPENDS test_material_standard_linear_isotropic_hardening_mesh
  FILES_TO_COPY test_material_elasto_plastic_linear_isotropic_hardening.dat
  PACKAGE core implicit
  )


#===============================================================================
