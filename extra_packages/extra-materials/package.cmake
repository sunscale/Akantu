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

package_declare(extra_materials
  DESCRIPTION "Add the extra list of materials in Akantu"
  DEPENDS lapack)

package_declare_sources(extra_materials
  material_extra_includes.hh

  material_damage/material_brittle.cc
  material_damage/material_brittle.hh
  material_damage/material_brittle_inline_impl.hh

  material_damage/material_damage_iterative.cc
  material_damage/material_damage_iterative.hh
  material_damage/material_damage_iterative_inline_impl.hh

  material_damage/material_iterative_stiffness_reduction.cc
  material_damage/material_iterative_stiffness_reduction.hh

  material_damage/material_damage_linear.cc
  material_damage/material_damage_linear.hh
  material_damage/material_damage_linear_inline_impl.hh

  material_damage/material_vreepeerlings.hh
  material_damage/material_vreepeerlings_inline_impl.hh
  material_damage/material_vreepeerlings_tmpl.hh

  material_plastic/material_viscoplastic.cc
  material_plastic/material_viscoplastic.hh
  material_plastic/material_viscoplastic_inline_impl.hh

  material_viscoelastic/material_stiffness_proportional.cc
  material_viscoelastic/material_stiffness_proportional.hh

  material_damage/material_orthotropic_damage.hh
  material_damage/material_orthotropic_damage_tmpl.hh

  material_damage/material_orthotropic_damage_iterative.cc
  material_damage/material_orthotropic_damage_iterative.hh
  material_damage/material_orthotropic_damage_iterative_inline_impl.hh

  material_FE2/material_FE2.hh
  material_FE2/material_FE2.cc
  material_FE2/material_FE2_inline_impl.hh
  material_FE2/solid_mechanics_model_RVE.hh
  material_FE2/solid_mechanics_model_RVE.cc
  )

package_declare_material_infos(extra_materials
  LIST AKANTU_EXTRA_MATERIAL_LIST
  INCLUDE material_extra_includes.hh
  )

package_declare_documentation_files(extra_materials
  manual-extra_materials.tex
  manual-appendix-materials-extra-materials.tex

  figures/stress_strain_visco.pdf
  )

package_declare_documentation(extra_materials
  "This package activates additional constitutive laws:"
  "\\begin{itemize}"
  "  \\item Linear anisotropy"
  "  \\item Linear orthotropy"
  "  \\item Visco-plastic"
  "\\end{itemize}"
  )

package_declare(extra_materials_non_local
  DESCRIPTION "Add the extra list of non local materials in Akantu"
  DEPENDS extra_materials damage_non_local)

package_declare_sources(extra_materials_non_local
  material_damage/material_orthotropic_damage_non_local.hh

  material_damage/material_vreepeerlings_non_local.cc
  material_damage/material_vreepeerlings_non_local.hh
  material_damage/material_brittle_non_local.hh
  material_damage/material_damage_iterative_non_local.hh
  material_damage/material_damage_iterative_non_local.cc
  material_damage/material_orthotropic_damage_iterative_non_local.hh

  material_damage/material_vreepeerlings_non_local_inline_impl.hh
  material_damage/material_brittle_non_local_inline_impl.hh
  material_damage/material_damage_iterative_non_local_inline_impl.hh
  material_damage/material_orthotropic_damage_iterative_non_local_inline_impl.hh

  material_non_local_extra_includes.hh
  )

package_declare_material_infos(extra_materials_non_local
  LIST AKANTU_DAMAGE_NON_LOCAL_MATERIAL_EXTRA_LIST
  INCLUDE material_extra_includes.hh
)
