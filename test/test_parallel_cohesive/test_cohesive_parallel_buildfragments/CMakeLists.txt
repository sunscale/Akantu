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
  DEPENDENCIES test_cohesive_parallel_buildfragments_mesh
  FILES_TO_COPY material.dat
  DIRECTORIES_TO_CREATE paraview)

add_dependencies(test_cohesive_parallel_buildfragments test_cohesive_parallel_buildfragments_mesh)
