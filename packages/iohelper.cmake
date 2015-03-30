#===============================================================================
# @file   90_iohelper.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Tue Nov 29 2011
# @date last modification: Tue Sep 02 2014
#
# @brief  package description for iohelper
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

package_declare(IOHelper EXTERNAL
  DESCRIPTION "Add IOHelper support in akantu"
  SYSTEM OFF third-party/cmake/iohelper.cmake
  DEFAULT ON)

package_add_third_party_script_variable(IOHelper
  IOHELPER_VERSION "1.1")
package_add_third_party_script_variable(IOHelper
  IOHELPER_GIT "https://git.epfl.ch/repo/iohelper.git")

package_declare_sources(IOHelper
  io/dumper/dumpable_iohelper.hh
  io/dumper/dumper_iohelper.hh
  io/dumper/dumper_iohelper.cc
  io/dumper/dumper_paraview.hh
  io/dumper/dumper_text.cc
  io/dumper/dumper_text.hh
  io/dumper/dumper_paraview.cc
  io/dumper/dumper_homogenizing_field.hh
  io/dumper/dumper_type_traits.hh
  io/dumper/dumper_compute.hh
  io/dumper/dumper_nodal_field.hh
  io/dumper/dumper_quadrature_points_field.hh
  io/dumper/dumper_variable.hh
  io/dumper/dumper_iterator_helper.hh
  io/dumper/dumper_connectivity_field.hh
  io/dumper/dumper_padding_helper.hh
  io/dumper/dumper_elemental_field.hh
  io/dumper/dumper_element_iterator.hh
  io/dumper/dumper_material_internal_field.hh
  io/dumper/dumper_generic_elemental_field.hh
  io/dumper/dumper_generic_elemental_field_tmpl.hh
  )

package_declare_documentation(IOHelper
  "This package activates the IOHelper facilities withing Akantu. This is mandatory if you want to be able to output Paraview files"
  "as well as any Dumper within Akantu."
  )