#===============================================================================
# @file   package.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
#
# @brief  package description for extra materials list
#
# @section LICENSE
#
# Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
