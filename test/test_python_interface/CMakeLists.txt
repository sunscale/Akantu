#===============================================================================
# @file   CMakeLists.txt
#
# @author Fabian Barras <fabian.barras@epfl.ch>
# @author Lucas Frerot <lucas.frerot@epfl.ch>
#
# @date creation: Fri Sep 03 2010
# @date last modification: Mon Feb 05 2018
#
# @brief  Python Interface tests
#
# @section LICENSE
#
# Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# Akantu is free  software: you can redistribute it and/or  modify it under the terms  of the  GNU Lesser  General Public  License as published by  the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more details.
#
# You should  have received  a copy  of the GNU  Lesser General  Public License along with Akantu. If not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================

add_mesh(mesh_python mesh_dcb_2d.geo ORDER 2 DIM 2)

register_test(test_multiple_init
  SCRIPT test_multiple_init.py
  PYTHON
  DEPENDS mesh_python
  PACKAGE python_interface
  FILES_TO_COPY input_test.dat
  )

register_test(test_mesh_interface
  SCRIPT test_mesh_interface.py
  PYTHON
  DEPENDS mesh_python
  PACKAGE python_interface
  )

register_test(test_boundary_condition_functors
  SCRIPT test_boundary_condition_functors.py
  PYTHON
  DEPENDS mesh_python
  PACKAGE python_interface
  FILES_TO_COPY input_test.dat)
