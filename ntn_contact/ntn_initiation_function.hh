/**
 * @file   ntn_initiation_function.hh
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Mon Sep  2 14:31:00 2013
 *
 * @brief  initiation ntn and ntrf friction
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
// simtools
#include "parameter_reader.hh"
#include "ntrf_contact.hh"
#include "ntn_base_friction.hh"

__BEGIN_SIMTOOLS__

using namespace akantu;

NTNBaseFriction * initializeNTNFriction(NTNBaseContact * contact, 
					ParameterReader & data);

__END_SIMTOOLS__
