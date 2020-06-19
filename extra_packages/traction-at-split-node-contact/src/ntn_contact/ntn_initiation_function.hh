/**
 * @file   ntn_initiation_function.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Fri Jan 04 2013
 * @date last modification: Fri Feb 23 2018
 *
 * @brief  initiation ntn and ntrf friction
 *
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
// simtools
#include "ntn_base_friction.hh"
#include "ntrf_contact.hh"
#include "parameter_reader.hh"

namespace akantu {

std::unique_ptr<NTNBaseFriction>
initializeNTNFriction(NTNBaseContact & contact);

std::unique_ptr<NTNBaseFriction>
initializeNTNFriction(NTNBaseContact & contact,
                      const std::string & friction_law,
                      const std::string & friction_reg);

} // namespace akantu
