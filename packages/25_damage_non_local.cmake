#===============================================================================
# @file   25_damage_non_local.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Fri Jun 15 2012
# @date last modification: Fri Jun 13 2014
#
# @brief  package description for non-local materials
#
# @section LICENSE
#
# Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

option(AKANTU_DAMAGE_NON_LOCAL "Package for Non-local damage constitutives laws Akantu" OFF)

add_external_package_dependencies(damage_non_local lapack)

set(AKANTU_DAMAGE_NON_LOCAL_FILES
  model/solid_mechanics/materials/material_damage/material_damage_non_local.hh
  model/solid_mechanics/materials/material_damage/material_marigo_non_local.hh
  model/solid_mechanics/materials/material_damage/material_marigo_non_local_inline_impl.cc
  model/solid_mechanics/materials/material_damage/material_mazars_non_local.cc
  model/solid_mechanics/materials/material_damage/material_mazars_non_local.hh

  model/solid_mechanics/materials/material_non_local.hh
  model/solid_mechanics/materials/material_non_local_includes.hh
  model/solid_mechanics/materials/material_non_local_inline_impl.cc

  model/solid_mechanics/materials/weight_function.cc
  model/solid_mechanics/materials/weight_function.hh
  model/solid_mechanics/materials/weight_function_tmpl.hh

  synchronizer/grid_synchronizer.cc
  synchronizer/grid_synchronizer.hh
  )

set(AKANTU_DAMAGE_NON_LOCAL_TESTS
  test_material_damage_non_local
  test_grid_synchronizer
  )

set(AKANTU_DAMAGE_NON_LOCAL_MANUAL_FILES
  manual-constitutive-laws-non_local.tex
  manual-appendix-materials-non-local.tex)

set(AKANTU_DAMAGE_NON_LOCAL_DOCUMENTATION "
This package activates the non local damage feature of AKANTU
")
