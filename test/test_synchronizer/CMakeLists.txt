#===============================================================================
# @file   CMakeLists.txt
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Fri Sep 03 2010
# @date last modification: Wed Jan 20 2016
#
# @brief  configuration for synchronizer tests
#
# @section LICENSE
#
# Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
# Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
# Solides)
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

add_mesh(test_synchronizer_communication_mesh
  cube.geo 3 2)

register_test(test_synchronizer_communication
  SOURCES test_synchronizer_communication.cc
  DEPENDS test_synchronizer_communication_mesh
  PACKAGE parallel
  )

# if(DEFINED AKANTU_DAMAGE_NON_LOCAL)
#   add_executable(test_grid_synchronizer_check_neighbors test_grid_synchronizer_check_neighbors.cc test_grid_tools.hh)
#   target_link_libraries(test_grid_synchronizer_check_neighbors akantu)
#   if(AKANTU_EXTRA_CXX_FLAGS)
#     set_target_properties(test_grid_synchronizer_check_neighbors PROPERTIES COMPILE_FLAGS ${AKANTU_EXTRA_CXX_FLAGS})
#   endif()
# endif()

# register_test(test_grid_synchronizer
#   SOURCES test_grid_synchronizer.cc test_data_accessor.hh
#   DEPENDS test_synchronizer_communication_mesh test_grid_synchronizer_check_neighbors
#   EXTRA_FILES test_grid_synchronizer_check_neighbors.cc test_grid_tools.hh
#   PACKAGE damage_non_local
#   )

# register_test(test_dof_synchronizer
#   SOURCES test_dof_synchronizer.cc test_data_accessor.hh
#   FILES_TO_COPY bar.msh
#   PACKAGE parallel
#   PARALLEL
#   )

# register_test(test_dof_synchronizer_communication
#   SOURCES test_dof_synchronizer_communication.cc test_dof_data_accessor.hh
#   DEPENDS test_synchronizer_communication_mesh
#   PACKAGE parallel
#   )

# register_test(test_data_distribution
#   SOURCES test_data_distribution.cc
#   FILES_TO_COPY data_split.msh
#   PACKAGE parallel
#   )
