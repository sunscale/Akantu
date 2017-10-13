#===============================================================================
# @file   CMakeLists.txt
#
# @author Seyedeh Mohadeseh Taheri Mousavi <mohadeseh.taherimousavi@epfl.ch>
# @author Marco Vocialta <marco.vocialta@epfl.ch>
#
#
# @brief  Configuration for test_cohesive_parallel_extrinsic_IG_TG
#
# @section LICENSE
#
# Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# @section DESCRIPTION
#
#===============================================================================

register_test(test_cohesive_parallel_extrinsic_IG_TG
  SOURCES test_cohesive_parallel_extrinsic_IG_TG.cc
  PACKAGE parallel_cohesive_element
  FILES_TO_COPY material.dat square.msh
  DIRECTORIES_TO_CREATE paraview)
