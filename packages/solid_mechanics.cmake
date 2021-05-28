#===============================================================================
# @file   solid_mechanics.cmake
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Mon Nov 21 2011
# @date last modification: Mon Jan 18 2016
#
# @brief  package description for core
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
#===============================================================================
package_declare(solid_mechanics DEFAULT ON
  DESCRIPTION "Solid mechanics model"
  DEPENDS core lapack
  )

package_declare_sources(solid_mechanics
  model/solid_mechanics/material.cc
  model/solid_mechanics/material.hh
  model/solid_mechanics/material_inline_impl.hh
  model/solid_mechanics/material_selector.hh
  model/solid_mechanics/material_selector_tmpl.hh
  model/solid_mechanics/materials/internal_field.hh
  model/solid_mechanics/materials/internal_field_tmpl.hh
  model/solid_mechanics/materials/random_internal_field.hh
  model/solid_mechanics/materials/random_internal_field_tmpl.hh
  model/solid_mechanics/solid_mechanics_model.cc
  model/solid_mechanics/solid_mechanics_model.hh
  model/solid_mechanics/solid_mechanics_model_inline_impl.hh
  model/solid_mechanics/solid_mechanics_model_io.cc
  model/solid_mechanics/solid_mechanics_model_mass.cc
  model/solid_mechanics/solid_mechanics_model_material.cc
  model/solid_mechanics/solid_mechanics_model_tmpl.hh
  model/solid_mechanics/solid_mechanics_model_event_handler.hh
  model/solid_mechanics/materials/plane_stress_toolbox.hh
  model/solid_mechanics/materials/plane_stress_toolbox_tmpl.hh

  model/solid_mechanics/materials/material_core_includes.hh
  model/solid_mechanics/materials/material_elastic.cc
  model/solid_mechanics/materials/material_elastic.hh
  model/solid_mechanics/materials/material_elastic_inline_impl.hh
  model/solid_mechanics/materials/material_thermal.cc
  model/solid_mechanics/materials/material_thermal.hh
  model/solid_mechanics/materials/material_elastic_linear_anisotropic.cc
  model/solid_mechanics/materials/material_elastic_linear_anisotropic.hh
  model/solid_mechanics/materials/material_elastic_linear_anisotropic_inline_impl.hh
  model/solid_mechanics/materials/material_elastic_orthotropic.cc
  model/solid_mechanics/materials/material_elastic_orthotropic.hh
  model/solid_mechanics/materials/material_damage/material_anisotropic_damage.hh
  model/solid_mechanics/materials/material_damage/material_anisotropic_damage.cc
  model/solid_mechanics/materials/material_damage/material_anisotropic_damage_tmpl.hh
  model/solid_mechanics/materials/material_damage/material_damage.hh
  model/solid_mechanics/materials/material_damage/material_damage_tmpl.hh
  model/solid_mechanics/materials/material_damage/material_marigo.cc
  model/solid_mechanics/materials/material_damage/material_marigo.hh
  model/solid_mechanics/materials/material_damage/material_marigo_inline_impl.hh
  model/solid_mechanics/materials/material_damage/material_mazars.cc
  model/solid_mechanics/materials/material_damage/material_mazars.hh
  model/solid_mechanics/materials/material_damage/material_phasefield.cc
  model/solid_mechanics/materials/material_damage/material_phasefield.hh
  model/solid_mechanics/materials/material_damage/material_phasefield_inline_impl.cc
  model/solid_mechanics/materials/material_damage/material_mazars_inline_impl.hh
  model/solid_mechanics/materials/material_finite_deformation/material_neohookean.cc
  model/solid_mechanics/materials/material_finite_deformation/material_neohookean.hh
  model/solid_mechanics/materials/material_finite_deformation/material_neohookean_inline_impl.hh
  model/solid_mechanics/materials/material_plastic/material_plastic.cc
  model/solid_mechanics/materials/material_plastic/material_plastic.hh
  model/solid_mechanics/materials/material_plastic/material_plastic_inline_impl.hh
  model/solid_mechanics/materials/material_plastic/material_linear_isotropic_hardening.cc
  model/solid_mechanics/materials/material_plastic/material_linear_isotropic_hardening.hh
  model/solid_mechanics/materials/material_plastic/material_linear_isotropic_hardening_inline_impl.hh
  model/solid_mechanics/materials/material_viscoelastic/material_standard_linear_solid_deviatoric.cc
  model/solid_mechanics/materials/material_viscoelastic/material_standard_linear_solid_deviatoric.hh
  model/solid_mechanics/materials/material_viscoelastic/material_viscoelastic_maxwell.cc
  model/solid_mechanics/materials/material_viscoelastic/material_viscoelastic_maxwell.hh

  model/solid_mechanics/materials/material_non_local.hh
  model/solid_mechanics/materials/material_non_local_tmpl.hh
  model/solid_mechanics/materials/material_non_local_includes.hh
  )

package_declare_material_infos(solid_mechanics
  LIST AKANTU_CORE_MATERIAL_LIST
  INCLUDE material_core_includes.hh
  )

package_declare_documentation_files(solid_mechanics
  manual-solidmechanicsmodel.tex
  manual-constitutive-laws.tex
  manual-lumping.tex
  manual-appendix-materials.tex

  figures/dynamic_analysis.png
  figures/explicit_dynamic.pdf
  figures/explicit_dynamic.svg
  figures/static.pdf
  figures/static.svg
  figures/hooke_law.pdf
  figures/implicit_dynamic.pdf
  figures/implicit_dynamic.svg
  figures/problemDomain.pdf_tex
  figures/problemDomain.pdf
  figures/static_analysis.png
  figures/stress_strain_el.pdf
  figures/tangent.pdf
  figures/tangent.svg

  figures/stress_strain_neo.pdf
  figures/visco_elastic_law.pdf
  figures/isotropic_hardening_plasticity.pdf
  figures/stress_strain_visco.pdf
  )

package_declare_extra_files_to_package(solid_mechanics
  SOURCES
    model/solid_mechanics/material_list.hh.in
  )
