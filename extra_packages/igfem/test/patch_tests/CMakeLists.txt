#===============================================================================
# @file   CMakeLists.txt
#
# @author Aurelia Cuba Ramos <aurelia.cubaramos@epfl.ch>
# @author Clement Roux-Langlois <clement.roux@epfl.ch>
#
#
# @brief  configuration for IGFEM tests
#
# @section LICENSE
#
# Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# @section DESCRIPTION
#
#===============================================================================

add_mesh(test_mesh plate.geo 2 1 OUTPUT plate.msh)

register_test(test_igfem_triangle_5
  SOURCES test_igfem_triangle_5.cc
  DEPENDS test_mesh
  FILES_TO_COPY material.dat
  DIRECTORIES_TO_CREATE paraview
  PACKAGE IGFEM
  )

register_test(test_igfem_triangle_4
  SOURCES test_igfem_triangle_4.cc
  FILES_TO_COPY material.dat test_mesh.msh
  DIRECTORIES_TO_CREATE paraview
  PACKAGE IGFEM
  )

register_test(test_interface_position
  SOURCES test_interface_position.cc
  FILES_TO_COPY material_test_interface_position.dat test_interface_position.msh
  DIRECTORIES_TO_CREATE paraview
  PACKAGE IGFEM
  )
