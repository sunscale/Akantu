#===============================================================================
# @file   CMakeLists.txt
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
#
# @brief  configuration for synchronizer tests
#
# @section LICENSE
#
# Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# @section DESCRIPTION
#
#===============================================================================

add_akantu_test(test_parallel_cohesive "parallel_cohesive_test")

add_mesh(test_facet_synchronizer_mesh
  facet.geo 3 2)

register_test(test_facet_synchronizer
  SOURCES test_facet_synchronizer.cc test_data_accessor.hh
  DEPENDS test_facet_synchronizer_mesh
  PACKAGE parallel_cohesive_element
  )
