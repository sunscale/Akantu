#===============================================================================
# @file   petsc.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
# @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
# @author Aurelia Cuba Ramos <aurelia.cubaramos@epfl.ch>
#
# @date   Mon Nov 21 18:19:15 2011
#
# @brief  package description for PETSc support
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
#===============================================================================

add_optional_external_package(PETSc "Add PETSc support in akantu" OFF ARGS COMPONENT ARGS CXX)
add_internal_package_dependencies(petsc parallel)

set(AKANTU_PETSC_FILES
  solver/petsc_matrix.hh
  solver/petsc_matrix.cc
  solver/petsc_matrix_inline_impl.cc
  solver/solver_petsc.hh
  solver/solver_petsc.cc
  )

set(AKANTU_PETSC_TESTS
  test_petsc_matrix_profile
  )
