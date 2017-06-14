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

register_test(test_transfer_internals
  SOURCES test_transfer_internals.cc
  FILES_TO_COPY material_damage.dat fine_mesh.msh 
  DIRECTORIES_TO_CREATE paraview
  PACKAGE IGFEM 
  )

register_test(test_eigenstrain
  SOURCES test_eigenstrain.cc
  FILES_TO_COPY material_damage.dat fine_mesh.msh 
  DIRECTORIES_TO_CREATE paraview
  PACKAGE IGFEM extra_materials
  )

register_test(test_igfem_element_orientation
  SOURCES test_igfem_element_orientation.cc
  FILES_TO_COPY material.dat test_igfem_element_orientation.msh
  DIRECTORIES_TO_CREATE paraview
  PACKAGE IGFEM
  )

register_test(test_igfem_element_orientation_parallel
  SOURCES test_igfem_element_orientation.cc
  FILES_TO_COPY material.dat test_igfem_element_orientation.msh
  DIRECTORIES_TO_CREATE paraview
  PACKAGE IGFEM
  )

register_test(test_ASR_damge
  SOURCES test_ASR_damage.cc
  FILES_TO_COPY material_ASR.dat one_inclusion.msh
  DIRECTORIES_TO_CREATE paraview
  PACKAGE IGFEM extra_materials
  )

register_test(test_material_igfem_iterative_stiffness_reduction
  SOURCES test_material_igfem_iterative_stiffness_reduction.cc
  FILES_TO_COPY material_stiffness_reduction.dat regular_mesh.msh igfem_mesh.msh
  DIRECTORIES_TO_CREATE paraview
  PACKAGE IGFEM implicit extra_materials
  )

register_test(test_material_igfem_iterative_stiffness_reduction_damage_step_transfer
  SOURCES test_material_igfem_iterative_stiffness_reduction_damage_step_transfer.cc
  FILES_TO_COPY material_stiffness_reduction_2.dat test_damage_transfer.msh
  DIRECTORIES_TO_CREATE paraview
  PACKAGE IGFEM implicit extra_materials
  )

register_test(test_material_igfem_iterative_stiffness_reduction_parallel
  SOURCES test_material_igfem_iterative_stiffness_reduction_parallel.cc
  FILES_TO_COPY material_stiffness_reduction_2.dat test_damage_transfer.msh
  DIRECTORIES_TO_CREATE paraview
  PACKAGE IGFEM implicit extra_materials
  )

register_test(test_material_igfem_iterative_stiffness_reduction_tangent_transfer
  SOURCES test_material_igfem_iterative_stiffness_reduction_tangent_transfer.cc
  FILES_TO_COPY material_stiffness_reduction_weibull.dat test_damage_transfer.msh
  DIRECTORIES_TO_CREATE paraview
  PACKAGE IGFEM implicit extra_materials
  )

register_test(test_volume_computation
  SOURCES test_volume_computation.cc
  FILES_TO_COPY material_stiffness_reduction_2.dat test_damage_transfer.msh
  DIRECTORIES_TO_CREATE paraview
  PACKAGE IGFEM implicit extra_materials
  )

# register_test(test_igfem_analysis
#   SOURCES igfem_analysis.cc
#   FILES_TO_COPY material_el_flipping.dat one_inclusion.msh 
#   DIRECTORIES_TO_CREATE paraview
#   PACKAGE IGFEM
#   )
# 
# register_test(test_element_flipping
#   SOURCES test_element_flipping.cc
#   FILES_TO_COPY standard_materials.dat
#   DIRECTORIES_TO_CREATE paraview
#   PACKAGE IGFEM
#   )
