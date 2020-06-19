/**
 * @file   material_cohesive_includes.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Sun Sep 26 2010
 * @date last modification: Fri Oct 13 2017
 *
 * @brief  List of includes for cohesive elements
 *
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

// /* --------------------------------------------------------------------------
// */
// #ifndef AKANTU_CMAKE_LIST_MATERIALS
// #include "material_cohesive.hh"
// #include "material_cohesive_bilinear.hh"
// #include "material_cohesive_exponential.hh"
// #include "material_cohesive_linear.hh"
// #include "material_cohesive_linear_fatigue.hh"
// #include "material_cohesive_linear_friction.hh"
// #include "material_cohesive_linear_uncoupled.hh"
// #endif

#define AKANTU_COHESIVE_MATERIAL_LIST                                          \
  ((2, (cohesive_linear, MaterialCohesiveLinear)))(                            \
      (2, (cohesive_linear_fatigue, MaterialCohesiveLinearFatigue)))(          \
      (2, (cohesive_linear_friction, MaterialCohesiveLinearFriction)))(        \
      (2, (cohesive_linear_uncoupled, MaterialCohesiveLinearUncoupled)))(      \
      (2, (cohesive_bilinear, MaterialCohesiveBilinear)))(                     \
      (2, (cohesive_exponential, MaterialCohesiveExponential)))
