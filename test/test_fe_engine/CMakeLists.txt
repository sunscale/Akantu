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

register_test(test_interpolate
  SOURCES test_interpolate.cc igfem_mesh_generation.cc
  FILES_TO_COPY mesh_info_igfem_triangle_5.txt mesh_info_igfem_triangle_4.txt
  PACKAGE IGFEM
  )

register_test(test_integrate
  SOURCES test_integrate.cc igfem_mesh_generation.cc
  FILES_TO_COPY mesh_info_igfem_triangle_5.txt mesh_info_igfem_triangle_4.txt
  PACKAGE IGFEM
  )

register_test(test_fe_engine_precomputation
  SOURCES test_fe_engine_precomputation.cc igfem_mesh_generation.cc
  FILES_TO_COPY mesh_info_igfem_triangle_5.txt mesh_info_igfem_triangle_4.txt
  PACKAGE IGFEM
  )

register_test(test_gradient
  SOURCES test_gradient.cc igfem_mesh_generation.cc
  FILES_TO_COPY mesh_info_igfem_triangle_5.txt mesh_info_igfem_triangle_4.txt
  PACKAGE IGFEM
  )
