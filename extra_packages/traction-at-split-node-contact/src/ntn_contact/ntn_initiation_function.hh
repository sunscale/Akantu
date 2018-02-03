/**
 * @file   ntn_initiation_function.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 *
 * @brief  initiation ntn and ntrf friction
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
// simtools
#include "parameter_reader.hh"
#include "ntrf_contact.hh"
#include "ntn_base_friction.hh"

namespace akantu {

NTNBaseFriction * initializeNTNFriction(NTNBaseContact * contact, 
					ParameterReader & data);

NTNBaseFriction * initializeNTNFriction(NTNBaseContact * contact);

NTNBaseFriction * initializeNTNFriction(NTNBaseContact * contact,
					const std::string & friction_law,
					const std::string & friction_reg);

} // namespace akantu
