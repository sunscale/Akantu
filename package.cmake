#===============================================================================
# @file   igfem.cmake
#
# @author Aurelia Cuba Ramos <aurelia.cubaramos@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date   Fri Mai 23 18:19:15 2011
#
# @brief  package description for interface-enriched generalized IGFEM
#
# @section LICENSE
#
# Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

option(AKANTU_IGFEM "Use Interface-enriched generalized FEM" OFF)

set(AKANTU_IGFEM_FILES
  fe_engine/element_class_igfem.hh
  fe_engine/shape_igfem.hh
  fe_engine/shape_igfem_inline_impl.cc
  fe_engine/element_class_igfem.cc
  fe_engine/element_classes/element_class_igfem_triangle_3_inline_impl.cc
  fe_engine/igfem_element.hh
  fe_engine/igfem_element.cc
  )

set(AKANTU_IGFEM_TESTS
  test_igfem_integrate
)

