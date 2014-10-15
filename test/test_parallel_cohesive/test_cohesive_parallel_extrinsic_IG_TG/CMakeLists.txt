#===============================================================================
# @file   CMakeLists.txt
#
# @author Seyedeh Mohadeseh Taheri Mousavi <mohadeseh.taherimousavi@epfl.ch>
# @author Marco Vocialta <marco.vocialta@epfl.ch>
#
# @date   Thu April 18 13:31:00 2013
#
# @brief  Configuration for test_cohesive_parallel_extrinsic_IG_TG
#
# @section LICENSE
#
# Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# Akantu is free  software: you can redistribute it and/or  modify it under the
# terms  of the  GNU Lesser  General Public  License as  published by  the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
# details.
#
# You should  have received  a copy  of the GNU  Lesser General  Public License
# along with Akantu. If not, see <http://www.gnu.org/licenses/>.
#
# @section DESCRIPTION
#
#===============================================================================

if(AKANTU_PARALLEL)
register_test(test_cohesive_parallel_extrinsic_IG_TG
  SOURCES test_cohesive_parallel_extrinsic_IG_TG.cc
  PACKAGE cohesive_parallel_element
  FILES_TO_COPY material.dat square.msh
  DIRECTORIES_TO_CREATE paraview)
endif()
