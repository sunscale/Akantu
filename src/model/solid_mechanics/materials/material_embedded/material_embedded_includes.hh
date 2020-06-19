/**
 * @file   material_embedded_includes.hh
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 *
 * @date creation: Fri Jan 04 2013
 * @date last modification: Fri Feb 09 2018
 *
 * @brief  List of includes for embedded elements
 *
 *
 * Copyright (©) 2014-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

/* -------------------------------------------------------------------------- */
#ifndef AKANTU_CMAKE_LIST_MATERIALS
#include "material_reinforcement.hh"
#endif

#define AKANTU_MATERIAL_REINFORCEMENT_LAW_TMPL_LIST                            \
  ((elastic, (MaterialElastic<1>)))(                                           \
      (plastic, (MaterialLinearIsotropicHardening<1>)))

#define AKANTU_EMBEDDED_MATERIAL_LIST                                          \
  ((2, (reinforcement, MaterialReinforcement)))
