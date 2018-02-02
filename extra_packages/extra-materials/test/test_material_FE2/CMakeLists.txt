#===============================================================================
# @file   CMakeLists.txt
#
# @author Aurelia Cuba Ramos <aurelia.cubaramos@epfl.ch>
#
#
# @brief  configuration for material_FE2 tests
#
# @section LICENSE
#
# Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# @section DESCRIPTION
#

#===============================================================================

register_test(test_material_FE2
  SOURCES test_material_FE2.cc
  FILES_TO_COPY material.dat periodic_plate.msh square.msh
  DIRECTORIES_TO_CREATE paraview
  PACKAGE extra_materials
  )

register_test(test_eigenstrain_homogenization
  SOURCES test_eigenstrain_homogenization.cc
  FILES_TO_COPY mesoscale_materials.dat one_inclusion.msh
  DIRECTORIES_TO_CREATE paraview
  PACKAGE extra_materials
  )

add_mesh(homogenized_plate homogenized_plate.geo 2 1)
register_test(test_elastic_homogenization
  SOURCES test_elastic_homogenization.cc
  DEPENDS homogenized_plate
  FILES_TO_COPY material_orthotropic.dat
  DIRECTORIES_TO_CREATE paraview
  PACKAGE implicit extra_materials
  )

add_mesh(periodic_plate periodic_plate.geo 2 1)
register_test(test_node_selection
  SOURCES  test_node_selection.cc 
  DEPENDS periodic_plate
  FILES_TO_COPY material_test_boundary.dat
  DIRECTORIES_TO_CREATE paraview
  PACKAGE implicit extra_materials
  )