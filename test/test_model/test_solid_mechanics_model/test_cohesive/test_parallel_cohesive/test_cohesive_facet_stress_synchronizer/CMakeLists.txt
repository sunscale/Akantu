#===============================================================================
# @file   CMakeLists.txt
#
# @author Marco Vocialta <marco.vocialta@epfl.ch>
#
#
# @brief  configuration for parallel test of facet stress synchronizer
#
# @section LICENSE
#
# Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# @section DESCRIPTION
#
#===============================================================================

add_mesh(test_cohesive_facet_stress_synchronizer_mesh tetrahedron.geo 3 2)

register_test(test_cohesive_facet_stress_synchronizer
  SOURCES test_cohesive_facet_stress_synchronizer.cc
  DEPENDS test_cohesive_facet_stress_synchronizer_mesh
  PACKAGE parallel_cohesive_element
  FILES_TO_COPY material.dat)
