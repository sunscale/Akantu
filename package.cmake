#===============================================================================
# @file   package.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
#
# @brief  package description for parallel cohesive elements
#
# @section LICENSE
#
# Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
#===============================================================================

package_declare(parallel_cohesive_element
  DESCRIPTION "Use parallel cohesive_element package of Akantu"
  DEPENDS cohesive_element parallel)

package_declare_sources(parallel_cohesive_element
  cohesive_element_inserter_parallel.cc
  cohesive_element_inserter_inline_impl.cc
  solid_mechanics_model_cohesive_parallel.hh
  solid_mechanics_model_cohesive_parallel.cc
  solid_mechanics_model_cohesive_parallel_inline_impl.cc
  facet_synchronizer.cc
  facet_synchronizer.hh
  facet_synchronizer_inline_impl.cc
  facet_stress_synchronizer.cc
  facet_stress_synchronizer.hh
  )

set(AKANTU_PARALLEL_COHESIVE_ELEMENT_TESTS
  test_cohesive_parallel_intrinsic
  test_cohesive_parallel_extrinsic
  test_cohesive_parallel_extrinsic_IG_TG
  test_cohesive_ghost_element_insertion
  test_facet_synchronizer
  test_cohesive_parallel_buildfragments
  )
