#===============================================================================
# @file   CMakeLists.txt
#
# @author Seyedeh Mohadeseh Taheri Mousavi <mohadeseh.taherimousavi@epfl.ch>
#
# @date creation: Mon Jan 18 2016
# @date last modification: Tue Jan 19 2016
#
# @brief  CMakeFile for new material example
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

#===============================================================================
add_mesh(new_local_material_barre_trou_mesh barre_trou.geo 2 2)

register_example(new_local_material
  SOURCES
    local_material_damage.cc
    local_material_damage.hh
    local_material_damage_inline_impl.hh
    new_local_material.cc
  DEPENDS
    new_local_material_barre_trou_mesh
  FILES_TO_COPY
    material.dat
  )

add_example(viscoelastic_maxwell  "Example on how to instantiate and obtaine energies from the viscoelastic maxwell material"      PACKAGE core)
