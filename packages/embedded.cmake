#===============================================================================
# @file   embedded.cmake
#
# @author Lucas Frérot <lucas.frerot@epfl.ch>
#
# @date creation: Tue Oct 16 2012
# @date last modification: Thu Jun 12 2014
#
# @brief  package descrition for embedded model use
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
#===============================================================================

package_declare(embedded 
  DESCRIPTION "Add support for the embedded solid mechanics model"
  DEPENDS CGAL)

package_declare_sources(embedded
  model/solid_mechanics/materials/material_embedded/material_embedded_includes.hh

  model/solid_mechanics/embedded_interface_intersector.hh
  model/solid_mechanics/embedded_interface_intersector.cc
  model/solid_mechanics/embedded_interface_model.hh
  model/solid_mechanics/embedded_interface_model.cc

  model/solid_mechanics/materials/material_embedded/embedded_internal_field.hh

  model/solid_mechanics/materials/material_embedded/material_reinforcement.hh
  model/solid_mechanics/materials/material_embedded/material_reinforcement.cc
  model/solid_mechanics/materials/material_embedded/material_reinforcement_inline_impl.cc

  model/solid_mechanics/materials/material_embedded/material_reinforcement_template.hh
  model/solid_mechanics/materials/material_embedded/material_reinforcement_template_inline_impl.cc
  )
