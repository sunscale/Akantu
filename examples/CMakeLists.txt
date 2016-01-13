#===============================================================================
# @file   CMakeLists.txt
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Fri Feb 24 2012
# @date last modification: Tue Sep 23 2014
#
# @brief  List of examples
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
# @section DESCRIPTION
#
#===============================================================================

include_directories(
  ${AKANTU_INCLUDE_DIRS}
  ${AKANTU_EXTERNAL_LIB_INCLUDE_DIR}
  )

#===============================================================================
add_example(new_material         "Example on how to add a new material in Akantu"      PACKAGE core)
add_example(boundary_conditions  "Example on hoy to set boundary conditions"           PACKAGE core)
add_example(explicit             "Example on how to run an explicit simulation"        PACKAGE core)
add_example(io                   "Example on how to perform Input/Output operations"   PACKAGE core)
add_example(implicit             "Example on how to run an implicit simulation"        PACKAGE implicit)
add_example(static               "Example on how to run a static simulation"           PACKAGE implicit)
add_example(parallel_2d          "Example of how to write a parallel code with Akantu" PACKAGE parallel)
add_example(cohesive_element     "Cohesive element examples"                           PACKAGE cohesive_element)
add_example(contact              "Examples on how to use contact within Akantu"        PACKAGE contact)
add_example(optimization         "Optimization examples"                               PACKAGE optimization)
add_example(structural_mechanics "Structural mechanics model examples"                 PACKAGE structural_mechanics)
add_example(heat_transfer        "Example on how to run heat transfer simulation"      PACKAGE heat_transfer)
add_example(embedded             "Example on how to run embedded model simulation"     PACKAGE embedded)
#===============================================================================
