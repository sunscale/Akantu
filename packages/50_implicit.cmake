#===============================================================================
# @file   50_implicit.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Tue Oct 16 2012
# @date last modification: Thu Jun 12 2014
#
# @brief  package description for the implicit solver
#
# @section LICENSE
#
# Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

add_meta_package(IMPLICIT "Add support for implicit time scheme" OFF SCOTCH MUMPS)

set(AKANTU_IMPLICIT_TESTS
  test_solid_mechanics_model_bar_traction2d_mass_not_lumped
  test_solid_mechanics_model_implicit_1d
  test_solid_mechanics_model_implicit_2d
  test_solid_mechanics_model_implicit_dynamic_2d
  )

set(AKANTU_IMPLICIT_DOCUMENTATION 
"
This package activates the sparse solver necessary to solve implicitely static/dynamic finite element problems.
It depends on:
\\begin{itemize}
\\item \\href{http://mumps.enseeiht.fr/}{MUMPS}, a parallel sparse direct solver.
\\item \\href{http://www.labri.fr/perso/pelegrin/scotch/}{Scotch}, a graph partitioner.
\\end{itemize}
")