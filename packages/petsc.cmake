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

package_declare(PETSc EXTERNAL
  DESCRIPTION "Add PETSc support in akantu"
  EXTRA_PACKAGE_OPTIONS ARGS COMPONENTS C
  DEPENDS parallel)

package_declare_sources(petsc
  model/dof_manager_petsc.hh
  model/dof_manager_petsc.cc
  solver/sparse_matrix_petsc.hh
  solver/sparse_matrix_petsc.cc
  solver/solver_petsc.hh
  solver/solver_petsc.cc
  solver/petsc_wrapper.hh
  )

