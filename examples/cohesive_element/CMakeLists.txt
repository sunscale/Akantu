#===============================================================================
# @file   CMakeLists.txt
#
# @author Seyedeh Mohadeseh Taheri Mousavi <mohadeseh.taherimousavi@epfl.ch>
#
# @date creation: Tue Jan 12 2016
# @date last modification: Mon Jan 18 2016
#
# @brief  CMakeLists for the cohesive examples
#
# @section LICENSE
#
# Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory
# (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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

add_example(cohesive_extrinsic "Extrinsic cohesive element" PACKAGE cohesive_element)
add_example(cohesive_intrinsic "Intrinsic cohesive element" PACKAGE cohesive_element )

add_example(cohesive_extrinsic_ig_tg
  "Extrinsic cohesive element with intergranular and transgranular material properties"
  PACKAGE cohesive_element)

#add_example(cohesive_extrinsic_implicit "Extrinsic cohesive element in implicit" PACKAGE cohesive_element)
