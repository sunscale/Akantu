#===============================================================================
# @file   implicit.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Tue Oct 16 2012
# @date last modification: Fri Aug 21 2015
#
# @brief  package description for the implicit solver
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
#===============================================================================
package_declare(implicit META
  DESCRIPTION "Add support for implicit time scheme")

set(AKANTU_IMPLICIT_SOLVER "Mumps"
  CACHE STRING "Solver activated in Akantu")
set_property(CACHE AKANTU_IMPLICIT_SOLVER PROPERTY STRINGS
  Mumps
  PETSc
  Mumps+PETSc
  )

if(AKANTU_IMPLICIT_SOLVER MATCHES "Mumps")
  package_add_dependencies(implicit PRIVATE Mumps)
else()
  package_remove_dependencies(implicit Mumps)
  set(AKANTU_USE_MUMPS OFF CACHE BOOL "" FORCE)
endif()

if(AKANTU_IMPLICIT_SOLVER MATCHES "PETSc")
  package_add_dependencies(implicit
    PRIVATE PETSc)
else()
  package_remove_dependency(implicit PETSc)
  set(AKANTU_USE_PETSC OFF CACHE BOOL "" FORCE)
endif()

package_declare_documentation(implicit
  "This package activates the sparse solver necessary to solve implicitely static/dynamic"
  "finite element problems."
  "It depends on:"
  "\\begin{itemize}"
  "  \\item \\href{http://mumps.enseeiht.fr/}{MUMPS}, a parallel sparse direct solver."
  "  \\item \\href{http://www.labri.fr/perso/pelegrin/scotch/}{Scotch}, a graph partitioner."
  "\\end{itemize}"
  )
