#===============================================================================
# @file   parallel_cohesive_element.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date   Tue Oct 16 14:05:02 2012
#
# @brief  package description for parallel cohesive elements
#
# @section LICENSE
#
# Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

option(AKANTU_PARALLEL_COHESIVE_ELEMENT "Use parallel cohesive_element package of Akantu" OFF)
add_internal_package_dependencies(parallel_cohesive_element cohesive_element parallel)

set(AKANTU_PARALLEL_COHESIVE_ELEMENT_FILES
  mesh_utils/cohesive_element_inserter_parallel.cc
  mesh_utils/cohesive_element_inserter_inline_impl.cc
  model/solid_mechanics/solid_mechanics_model_cohesive_parallel.hh
  model/solid_mechanics/solid_mechanics_model_cohesive_parallel.cc
  model/solid_mechanics/solid_mechanics_model_cohesive_inline_impl.cc
  synchronizer/facet_synchronizer.cc
  synchronizer/facet_synchronizer.hh
  synchronizer/facet_synchronizer_inline_impl.cc
  synchronizer/facet_stress_synchronizer.cc
  synchronizer/facet_stress_synchronizer.hh
  )

set(AKANTU_PARALLEL_COHESIVE_ELEMENT_TESTS
  test_cohesive_parallel_intrinsic
  test_cohesive_parallel_extrinsic
  test_cohesive_parallel_extrinsic_IG_TG
  test_cohesive_ghost_element_insertion
  test_facet_synchronizer
  test_cohesive_parallel_buildfragments
  )
