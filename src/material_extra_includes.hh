/**
 * @file   material_extra_includes.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Wed Oct 31 16:24:42 2012
 *
 * @brief  Extra list of materials
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#ifndef AKANTU_CMAKE_LIST_MATERIALS

// visco-elastic materials
#include "material_stiffness_proportional.hh"

// damage materials
#include "material_brittle.hh"
#include "material_damage_iterative.hh"
#include "material_damage_linear.hh"
#include "material_vreepeerlings.hh"

// plasticity
#include "material_viscoplastic.hh"

#endif

#define  AKANTU_EXTRA_MATERIAL_LIST                                     \
  ((2, (damage_linear      , MaterialDamageLinear                 )))   \
  ((2, (brittle            , MaterialBrittle                      )))   \
  ((2, (damage_iterative   , MaterialDamageIterative              )))   \
  ((2, (vreepeerlings      , MaterialVreePeerlings                )))   \
  ((2, (ve_stiffness_prop  , MaterialStiffnessProportional        )))   \
  ((2, (visco_plastic      , MaterialViscoPlastic                 )))
