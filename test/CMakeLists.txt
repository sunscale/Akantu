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

add_akantu_test(test_solid_mechanics_model_igfem "Test solid mechanics model igfem")
add_akantu_test(test_geometry "Test the IGFEM intersection algorithms")
add_akantu_test(test_fe_engine "Test the IGFEM interpolation and integration")
add_akantu_test(patch_tests "Patch tests for IGFEM model")