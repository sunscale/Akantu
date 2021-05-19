#===============================================================================
# @file   CMakeLists.txt
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Fri Oct 22 2010
# @date last modification: Fri Jan 22 2016
#
# @brief  List of examples
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

add_example(new_material         "Example on how to add a new material in Akantu"      PACKAGE core)
add_example(boundary_conditions  "Example on hoy to set boundary conditions"           PACKAGE core)
add_example(explicit             "Example on how to run an explicit simulation"        PACKAGE core)
add_example(io                   "Example on how to perform Input/Output operations"   PACKAGE core)
add_example(implicit             "Example on how to run an implicit simulation"        PACKAGE implicit)
add_example(static               "Example on how to run a static simulation"           PACKAGE implicit)
add_example(parallel             "Example of how to write a parallel code with Akantu" PACKAGE parallel)
add_example(cohesive_element     "Cohesive element examples"                           PACKAGE cohesive_element)
add_example(structural_mechanics "Structural mechanics model examples"                 PACKAGE structural_mechanics)
add_example(heat_transfer        "Example on how to run heat transfer simulation"      PACKAGE heat_transfer)
add_example(python               "Example on how to use the python interface"          PACKAGE python_interface)
add_example(embedded             "Example on how to run embedded model simulation"     PACKAGE embedded)
add_example(phase_field          "Example on how to run phase field model simulation"  PACKAGE phase_field)


package_add_files_to_package(
  examples/README.rst
  cmake/AkantuExampleMacros.cmake
  )


