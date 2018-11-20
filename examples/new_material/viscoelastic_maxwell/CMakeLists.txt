#===============================================================================
# @file   CMakeLists.txt
#
# @author Emil Gallyamov <emil.gallyamov@epfl.ch>
#
# @date creation: Tue Nov 20 2018
#
# @brief  CMakeFile for viscoelastic material example
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




add_mesh(material_viscoelastic_maxwell_mesh material_viscoelastic_maxwell_mesh.geo 2 1)
register_example(material_viscoelastic_maxwell_energies
  SOURCES material_viscoelastic_maxwell_energies.cc
  DEPENDS material_viscoelastic_maxwell_mesh
  USE_PACKAGES iohelper
  FILES_TO_COPY material_viscoelastic_maxwell.dat
  DIRECTORIES_TO_CREATE paraview
  )

