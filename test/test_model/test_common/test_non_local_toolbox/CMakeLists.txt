#===============================================================================
# @file   CMakeLists.txt
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Fri Sep 03 2010
# @date last modification: Fri Feb 02 2018
#
# @brief  configuration for heat transfer model tests
#
# @section LICENSE
#
# Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# Akantu is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with Akantu. If not, see <http://www.gnu.org/licenses/>.
#
# @section DESCRIPTION
#
#===============================================================================

#===============================================================================
package_is_activated(damage_non_local _act)
if(_act)
  add_library(test_material test_material.hh test_material.cc)
  add_library(test_material_damage test_material_damage.hh test_material_damage.cc)
  target_link_libraries(test_material akantu)
  target_link_libraries(test_material_damage akantu)
endif()

register_test(test_non_local_neighborhood_base
  SOURCES test_non_local_neighborhood_base.cc
  FILES_TO_COPY material.dat plot_neighborhoods.py plate.msh
  PACKAGE damage_non_local
  )

register_test(test_weight_computation
  SOURCES test_weight_computation.cc
  FILES_TO_COPY material_weight_computation.dat plate.msh
  PACKAGE damage_non_local
  )

register_test(test_non_local_averaging
  SOURCES test_non_local_averaging.cc
  FILES_TO_COPY material_avg.dat plate.msh
  PACKAGE damage_non_local
  LINK_LIBRARIES test_material
  )

register_test(test_remove_damage_weight_function
  SOURCES test_remove_damage_weight_function.cc
  FILES_TO_COPY material_remove_damage.dat plate.msh
  PACKAGE damage_non_local
  LINK_LIBRARIES test_material_damage
  )

register_test(test_build_neighborhood_parallel
  SOURCES test_build_neighborhood_parallel.cc
  FILES_TO_COPY material_parallel_test.dat parallel_test.msh
  PARALLEL
  PARALLEL_LEVEL 2
  PACKAGE damage_non_local
  LINK_LIBRARIES test_material
  )

register_test(test_pair_computation
  SOURCES test_pair_computation.cc
  FILES_TO_COPY material_remove_damage.dat pair_test.msh
  PARALLEL
  PARALLEL_LEVEL 1 2
  PACKAGE damage_non_local
  LINK_LIBRARIES test_material_damage
  )
