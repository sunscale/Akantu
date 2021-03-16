#===============================================================================
# @file   CMakeLists.txt
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Fri Sep 03 2010
# @date last modification: Mon Feb 12 2018
#
# @brief  configuration for tests
#
# @section LICENSE
#
# Copyright (©) 2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# Akantu is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with Akantu. If not, see <http://www.gnu.org/licenses/>.
#
# @section DESCRIPTION
#
#===============================================================================

include_directories(
  ${AKANTU_INCLUDE_DIRS}
  ${AKANTU_EXTERNAL_LIB_INCLUDE_DIR}
  )

set(AKANTU_TESTS_FILES CACHE INTERNAL "")

#===============================================================================
# List of tests
#===============================================================================
add_akantu_test(test_common "Test the common part of Akantu")
add_akantu_test(test_fe_engine "Test finite element functionalties")
add_akantu_test(test_mesh_utils "Test mesh utils")
add_akantu_test(test_mesh "Test mesh")
add_akantu_test(test_model "Test model objects")
add_akantu_test(test_solver "Test solver function")
add_akantu_test(test_io "Test the IO modules")
add_akantu_test(test_contact "Test the contact part of Akantu")
add_akantu_test(test_geometry "Test the geometry module of Akantu")
add_akantu_test(test_synchronizer "Test synchronizers")
add_akantu_test(test_python_interface "Test python interface")


package_add_files_to_package(
  cmake/akantu_test_driver.sh
  cmake/AkantuTestsMacros.cmake
  )

package_is_activated(parallel _is_parallel)
if (_is_parallel)
  option(AKANTU_TESTS_ALWAYS_USE_MPI "Defines if sequential tests should also use MPIEXEC" FALSE)
  mark_as_advanced(AKANTU_TESTS_ALWAYS_USE_MPI)
endif()

package_is_activated(gbenchmark _has_gbenchmark)
if (_has_gbenchmark)
  add_subdirectory(benchmark)
endif()
