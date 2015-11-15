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

register_test(test_spherical_kernel
  SOURCES test_spherical_kernel.cc
  PACKAGE CGAL
  )

register_test(test_generation_igfem_mesh_from_sphere_static
  SOURCES test_generation_igfem_mesh_from_sphere_static.cc
  FILES_TO_COPY test_geometry_triangle.msh
  PACKAGE CGAL IGFEM
  )

register_test(test_generation_igfem_mesh_from_sphere
  SOURCES test_generation_igfem_mesh_from_sphere.cc
  FILES_TO_COPY test_geometry_triangle.msh
  PACKAGE CGAL IGFEM
  )
