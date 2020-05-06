#===============================================================================
# @file   damage_non_local.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Fri Jun 15 2012
# @date last modification: Mon Jan 18 2016
#
# @brief  package description for non-local materials
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

package_declare(damage_non_local
  DESCRIPTION "Package for Non-local damage constitutives laws Akantu"
  DEPENDS lapack)

package_declare_sources(damage_non_local
  model/solid_mechanics/materials/material_damage/material_damage_non_local.hh
  model/solid_mechanics/materials/material_damage/material_marigo_non_local.cc
  model/solid_mechanics/materials/material_damage/material_marigo_non_local.hh
  model/solid_mechanics/materials/material_damage/material_mazars_non_local.cc
  model/solid_mechanics/materials/material_damage/material_mazars_non_local.hh

  model/solid_mechanics/materials/weight_functions/damaged_weight_function.hh
  model/solid_mechanics/materials/weight_functions/damaged_weight_function_inline_impl.hh
  model/solid_mechanics/materials/weight_functions/remove_damaged_weight_function.hh
  model/solid_mechanics/materials/weight_functions/remove_damaged_weight_function_inline_impl.hh
  model/solid_mechanics/materials/weight_functions/remove_damaged_with_damage_rate_weight_function.hh
  model/solid_mechanics/materials/weight_functions/remove_damaged_with_damage_rate_weight_function_inline_impl.hh
  model/solid_mechanics/materials/weight_functions/stress_based_weight_function.hh
  model/solid_mechanics/materials/weight_functions/stress_based_weight_function.cc
  model/solid_mechanics/materials/weight_functions/stress_based_weight_function_inline_impl.hh
  )

package_declare_material_infos(damage_non_local
  LIST AKANTU_DAMAGE_NON_LOCAL_MATERIAL_LIST
  INCLUDE material_non_local_includes.hh
  )

package_declare_documentation_files(damage_non_local
  manual-constitutive-laws-non_local.tex
  manual-appendix-materials-non-local.tex)

package_declare_documentation(damage_non_local
"This package activates the non local damage feature of AKANTU"
"")

