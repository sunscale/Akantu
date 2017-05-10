#===============================================================================
# @file   parallel.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Tue Oct 16 2012
# @date last modification: Fri Jul 10 2015
#
# @brief  meta package description for parallelization
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

package_declare(parallel META
  DESCRIPTION "Add parallel support in Akantu"
  DEPENDS mpi scotch)

set(AKANTU_PARALLEL_TESTS
  test_solid_mechanics_model_bar_traction2d_parallel
  test_solid_mechanics_model_segment_parallel
  test_solid_mechanics_model_pbc_parallel
  test_synchronizer_communication
  test_node_synchronizer
  test_dof_synchronizer
  test_dof_synchronizer_communication
  test_solid_mechanics_model_reassign_material
  )

package_declare_documentation_files(parallel
  manual-parallel.tex
  )

package_declare_documentation(parallel
  "This option activates the parallel features of AKANTU.")
