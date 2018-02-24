#===============================================================================
# @file   CMakeLists.txt
#
# @author Fabian Barras <fabian.barras@epfl.ch>
# @author Mauro Corrado <mauro.corrado@epfl.ch>
# @author Marco Vocialta <marco.vocialta@epfl.ch>
#
# @date creation: Wed Nov 05 2014
# @date last modification: Fri Oct 13 2017
#
# @brief  Tests parallel insertion of cohesive elements
#
# @section LICENSE
#
# Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# Akantu is free  software: you can redistribute it and/or  modify it under the terms  of the  GNU Lesser  General Public  License as published by  the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more details.
#
# You should  have received  a copy  of the GNU  Lesser General  Public License along with Akantu. If not, see <http://www.gnu.org/licenses/>.
#
# @section DESCRIPTION
#
#===============================================================================

add_mesh(3d_spherical_inclusion_parallel 3d_spherical_inclusion.geo 3 2)

register_test(test_cohesive_parallel_insertion_along_physical_surfaces
  SOURCES test_cohesive_parallel_insertion_along_physical_surfaces.cc
  DEPENDS 3d_spherical_inclusion_parallel
  FILES_TO_COPY input_file.dat
  PACKAGE parallel_cohesive_element
  )

register_test(test_cohesive_parallel_intrinsic_implicit_insertion
  PARALLEL
  POSTPROCESS verify_insertion.sh 
  PARALLEL_LEVEL "2 4"
  SOURCES test_cohesive_parallel_intrinsic_implicit_insertion.cc
  FILES_TO_COPY input_file_iii.dat
  2d_basic_interface.msh
  output_dir_verified
  PACKAGE parallel_cohesive_element
  DIRECTORIES_TO_CREATE output_dir  
  )

