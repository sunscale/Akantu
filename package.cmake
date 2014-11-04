#===============================================================================
# @file   extra_materials.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date   Wed Oct 31 16:24:42 2012
#
# @brief  package description for extra materials list
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

package_declare(extra-materials
  DESCRIPTION "Add the extra list of materials in Akantu"
  DEPENDS lapack)

package_declare_sources(extra_materials
  material_extra_includes.hh

  material_damage/material_brittle.cc
  material_damage/material_brittle.hh
  material_damage/material_brittle_inline_impl.cc

  material_damage/material_damage_iterative.cc
  material_damage/material_damage_iterative.hh
  material_damage/material_damage_iterative_inline_impl.cc

  material_damage/material_damage_linear.cc
  material_damage/material_damage_linear.hh
  material_damage/material_damage_linear_inline_impl.cc

  material_damage/material_vreepeerlings.hh
  material_damage/material_vreepeerlings_inline_impl.cc
  material_damage/material_vreepeerlings_tmpl.hh

  material_plastic/material_viscoplastic.cc
  material_plastic/material_viscoplastic.hh
  material_plastic/material_viscoplastic_inline_impl.cc

  material_viscoelastic/material_stiffness_proportional.cc
  material_viscoelastic/material_stiffness_proportional.hh

  material_damage/material_vreepeerlings_non_local.cc
  material_damage/material_vreepeerlings_non_local.hh
  material_damage/material_brittle_non_local.hh
  material_damage/material_damage_iterative_non_local.hh

  material_damage/material_vreepeerlings_non_local_inline_impl.cc
  material_damage/material_brittle_non_local_inline_impl.cc
  material_damage/material_damage_iterative_non_local_inline_impl.cc

  material_non_local_extra_includes.hh
  )

set(AKANTU_EXTRA_MATERIALS_TESTS
  )

set(AKANTU_EXTRA_MATERIALS_MANUAL_FILES
  manual-extra_materials.tex
  manual-appendix-materials-extra-materials.tex

  figures/stress_strain_visco.pdf
  )

set(AKANTU_EXTRA_MATERIALS_DOCUMENTATION "
This package activates additional constitutive laws:
\\begin{itemize}
\\item Linear anisotropy
\\item Linear orthotropy
\\item Visco-plastic
\\end{itemize}
" )
