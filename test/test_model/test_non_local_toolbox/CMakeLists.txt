#===============================================================================
# @file   CMakeLists.txt
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
#
# @date creation: Sat Sep 26 2015
# @date last modification: Wed Dec 16 2015
#
# @brief  configuration for heat transfer model tests
#
# @section LICENSE
#
# Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory
# (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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

#===============================================================================
register_test(test_non_local_neighborhood_base
  SOURCES test_non_local_neighborhood_base.cc test_material.hh test_material.cc
  FILES_TO_COPY material.dat plot_neighborhoods.py plate.msh
  PACKAGE damage_non_local
  )

register_test(test_weight_computation
  SOURCES test_weight_computation.cc test_material.hh test_material.cc 
  FILES_TO_COPY material_weight_computation.dat plate.msh
  PACKAGE damage_non_local
  )

register_test(test_non_local_averaging
  SOURCES test_non_local_averaging.cc test_material.hh test_material.cc 
  FILES_TO_COPY material_weight_computation.dat plate.msh
  PACKAGE damage_non_local
  )

register_test(test_remove_damage_weight_function
  SOURCES test_remove_damage_weight_function.cc test_material_damage.hh test_material_damage.cc 
  FILES_TO_COPY material_remove_damage.dat plate.msh
  PACKAGE damage_non_local
  )

register_test(test_build_neighborhood_parallel
  SOURCES test_build_neighborhood_parallel.cc test_material.hh test_material.cc
  FILES_TO_COPY material_parallel_test.dat parallel_test.msh
  PARALLEL
  PARALLEL_LEVEL 2
  PACKAGE damage_non_local
  )

register_test(test_pair_computation
  SOURCES test_pair_computation.cc test_material_damage.hh test_material_damage.cc
  FILES_TO_COPY material_remove_damage.dat pair_test.msh
  PARALLEL
  PARALLEL_LEVEL 1 2
  PACKAGE damage_non_local
  )
