#===============================================================================
# @file   CMakeLists.txt
#
# @author Marco Vocialta <marco.vocialta@epfl.ch>
#
#
# @brief  configuration for build fragments tests
#
# @section LICENSE
#
# Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# @section DESCRIPTION
#
#===============================================================================

add_mesh(test_cohesive_parallel_buildfragments_mesh mesh.geo 3 2)

register_test(test_cohesive_parallel_buildfragments
  SOURCES test_cohesive_parallel_buildfragments.cc
  DEPENDS test_cohesive_parallel_buildfragments_mesh
  FILES_TO_COPY material.dat
  PACKAGE parallel_cohesive_element
  DIRECTORIES_TO_CREATE paraview)
