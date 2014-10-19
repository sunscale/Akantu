#===============================================================================
# @file   CMakeLists.txt
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
# @author Cyprien Wolff <cyprien.wolff@epfl.ch>
#
# @date creation: Mon Apr 16 2012
# @date last modification: Thu Nov 01 2012
#
# @brief  configuration for explicit example
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

register_example(explicit_dynamic explicit_dynamic.cc)

add_mesh(explicit_dynamic_mesh bar.geo 3 1)
add_dependencies(explicit_dynamic explicit_dynamic_mesh)

#===============================================================================
file(COPY material.dat DESTINATION .)
