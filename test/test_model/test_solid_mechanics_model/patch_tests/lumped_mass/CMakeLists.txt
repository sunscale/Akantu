#===============================================================================
# @file   CMakeLists.txt
#
# @author David Simon Kammer <david.kammer@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Thu Mar 27 2014
# @date last modification: Thu Mar 27 2014
#
# @brief  configuration for the explicit patch test
#
# @section LICENSE
#
# Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

macro(register_patch_test name type)
  register_test(test_lumped_mass_${name}
    SOURCES test_lumped_mass.cc
    FILES_TO_COPY ../data/_${type}.msh ../data/material.dat
    COMPILE_OPTIONS "TYPE=_${type}"
    PACKAGE core
    )
endmacro()

foreach(_type ${LIST_TYPES} ${LIST_TYPES_2D})
  register_patch_test(${_type} ${_type})
endforeach()
