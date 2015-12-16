#===============================================================================
# @file   CMakeLists.txt
#
# @author Aurelia Cuba Ramos <aurelia.cubaramos@epfl.ch>
#
#
# @brief  configuration for material_damage tests
#
# @section LICENSE
#
# Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# @section DESCRIPTION
#

#===============================================================================

register_test(test_material_damage_iterative
  SOURCES test_material_damage_iterative.cc
  FILES_TO_COPY material.dat plate.msh
  DIRECTORIES_TO_CREATE paraview
  PACKAGE extra_materials
  )

register_test(test_material_damage_iterative_parallel
  SOURCES test_material_damage_iterative.cc
  FILES_TO_COPY material.dat plate.msh
  DIRECTORIES_TO_CREATE paraview
  PACKAGE extra_materials
  )

register_test(test_material_damage_iterative_non_local_serial
  SOURCES test_material_damage_iterative_non_local_serial.cc
  FILES_TO_COPY material_non_local.dat plate.msh
  DIRECTORIES_TO_CREATE paraview
  PACKAGE extra_materials_non_local
  )

register_test(test_material_damage_iterative_non_local_parallel
  SOURCES test_material_damage_iterative_non_local_parallel.cc
  FILES_TO_COPY two_materials.dat one_circular_inclusion.msh
  DIRECTORIES_TO_CREATE paraview
  PACKAGE extra_materials_non_local
  )

