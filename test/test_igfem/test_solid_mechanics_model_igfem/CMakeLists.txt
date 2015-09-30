#===============================================================================
# @file   CMakeLists.txt
#
# @author Aurelia Cuba Ramos <aurelia.cubaramos@epfl.ch>
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

register_test(test_solid_mechanics_model_igfem
  SOURCES test_solid_mechanics_model_igfem.cc
  FILES_TO_COPY material.dat plate.msh
  DIRECTORIES_TO_CREATE paraview
  PACKAGE IGFEM implicit
  )

register_test(test_moving_interface
  SOURCES test_moving_interface.cc
  FILES_TO_COPY material.dat plate.msh
  DIRECTORIES_TO_CREATE paraview
  PACKAGE IGFEM implicit
  )

register_test(test_transfer_internals
  SOURCES test_internal_transfer.cc
  FILES_TO_COPY material_damage.dat 
  DIRECTORIES_TO_CREATE paraview
  PACKAGE IGFEM 
  )

register_test(test_ASR_damge
  SOURCES test_ASR_damage.cc
  FILES_TO_COPY material_ASR.dat 
  DIRECTORIES_TO_CREATE paraview
  PACKAGE IGFEM extra_materials
  )
