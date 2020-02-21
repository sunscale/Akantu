#===============================================================================
# @file   CMakeLists.txt
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Fri Sep 03 2010
# @date last modification: Tue Dec 05 2017
#
# @brief  configurations for common tests
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
#===============================================================================

add_mesh(test_grid_mesh circle.geo 2 1)

register_test(test_grid test_grid.cc
  DEPENDS test_grid_mesh
  PACKAGE core)

#register_test(test_types test_types.cc PACKAGE core)

register_gtest_sources(SOURCES test_csr.cc PACKAGE core)
register_gtest_sources(SOURCES test_iterators.cc PACKAGE core HEADER_ONLY)
register_gtest_sources(SOURCES test_array.cc PACKAGE core)
register_gtest_sources(SOURCES test_tensors.cc PACKAGE core)
register_gtest_sources(SOURCES test_voigt_helper.cc PACKAGE core)

register_gtest_test(test_common)
