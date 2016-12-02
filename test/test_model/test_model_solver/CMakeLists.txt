#===============================================================================
# @file   CMakeLists.txt
#
# @brief  test for the common solvers interface of the models
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

register_test(test_dof_manager_default
  SOURCES test_dof_manager_default.cc
  PACKAGE mumps
  )

register_test(test_model_solver
  SOURCES test_model_solver.cc
  SCRIPT  test_model_solver.py
  PACKAGE mumps
  )

register_test(test_model_solver_dynamic_explicit
  SOURCES test_model_solver_dynamic.cc
  PACKAGE core
  COMPILE_OPTIONS "EXPLICIT=true"
  )

register_test(test_model_solver_dynamic_implicit
  SOURCES test_model_solver_dynamic.cc
  PACKAGE mumps
  COMPILE_OPTIONS "EXPLICIT=false"
  )
