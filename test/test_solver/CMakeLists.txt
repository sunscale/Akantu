#===============================================================================
# @file   CMakeLists.txt
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Fri Sep 03 2010
# @date last modification: Wed Dec 16 2015
#
# @brief  configuration for solver tests
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

add_mesh(test_solver_mesh triangle.geo 2 1)
add_mesh(test_matrix_mesh square.geo 2 1)
add_mesh(test_solver_petsc_mesh 1D_bar.geo 1 1)

register_test(test_sparse_matrix_profile
  SOURCES test_sparse_matrix_profile.cc
  DEPENDS test_solver_mesh
  PACKAGE implicit
  )

register_test(test_sparse_matrix_assemble
  SOURCES test_sparse_matrix_assemble.cc
  DEPENDS test_solver_mesh
  PACKAGE implicit
  )

register_test(test_sparse_matrix_product
  SOURCES test_sparse_matrix_product.cc
  FILES_TO_COPY bar.msh
  PACKAGE implicit
  )

register_test(test_sparse_solver_mumps
  SOURCES test_sparse_solver_mumps.cc
  PACKAGE mumps
  )

register_test(test_petsc_matrix_profile
  SOURCES test_petsc_matrix_profile.cc
  DEPENDS test_matrix_mesh
  PACKAGE petsc
  )

register_test(test_petsc_matrix_profile_parallel
  SOURCES test_petsc_matrix_profile_parallel.cc
  DEPENDS test_matrix_mesh
  PACKAGE petsc
  )

register_test(test_petsc_matrix_diagonal
  SOURCES test_petsc_matrix_diagonal.cc
  DEPENDS test_solver_mesh
  PACKAGE petsc
  )

register_test(test_petsc_matrix_apply_boundary
  SOURCES test_petsc_matrix_apply_boundary.cc
  DEPENDS test_solver_mesh
  PACKAGE petsc
  )

register_test(test_solver_petsc
  SOURCES test_solver_petsc.cc
  DEPENDS test_solver_petsc_mesh
  PACKAGE petsc
  )

register_test(test_solver_petsc_parallel
  SOURCES test_solver_petsc.cc
  DEPENDS test_solver_petsc_mesh
  PACKAGE petsc
  )
