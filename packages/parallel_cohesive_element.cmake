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
  mesh_utils/cohesive_element_inserter_parallel.cc
  mesh_utils/cohesive_element_inserter_inline_impl.cc
  model/solid_mechanics/solid_mechanics_model_cohesive/solid_mechanics_model_cohesive_parallel.hh
  model/solid_mechanics/solid_mechanics_model_cohesive/solid_mechanics_model_cohesive_parallel.cc
  model/solid_mechanics/solid_mechanics_model_cohesive/solid_mechanics_model_cohesive_parallel_inline_impl.cc
  )

package_declare_documentation(parallel_cohesive_element
  "This option activates the parallel cohesive elements' features of AKANTU.")
