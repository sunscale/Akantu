#===============================================================================
# @file   CMakeLists.txt
#
# @author Marco Vocialta <marco.vocialta@epfl.ch>
#
#
# @brief  configuration for cohesive elements tests
#
# @section LICENSE
#
# Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# @section DESCRIPTION
#
#===============================================================================

add_akantu_test(test_cohesive_parallel_intrinsic "test_cohesive_parallel_intrinsic")
add_akantu_test(test_cohesive_facet_stress_synchronizer "test_cohesive_facet_stress_synchronizer")
add_akantu_test(test_cohesive_parallel_extrinsic "test_cohesive_parallel_extrinsic")
add_akantu_test(test_cohesive_parallel_extrinsic_IG_TG "test_cohesive_parallel_extrinsic_IG_TG")
add_akantu_test(test_cohesive_parallel_buildfragments "test_cohesive_parallel_buildfragments")
add_akantu_test(test_cohesive_parallel_insertion "test_cohesive_parallel_insertion")
