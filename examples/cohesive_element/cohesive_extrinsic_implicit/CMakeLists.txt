#===============================================================================
# @file   CMakeLists.txt
#
# @author Mauro Corrado <mauro.corrado@epfl.ch>
#
# @date creation: Tue Jan 12 2016
#
# @brief  Example for extrinsic cohesive elements in implicit
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
# @section DESCRIPTION
#
#===============================================================================

add_mesh(cohesive_extrinsic_implicit_mesh dcb_2d.geo 2 2)

register_example(cohesive_extrinsic_implicit
  SOURCES cohesive_extrinsic_implicit.cc
  DEPENDS cohesive_extrinsic_implicit_mesh
  FILES_TO_COPY material.dat)

