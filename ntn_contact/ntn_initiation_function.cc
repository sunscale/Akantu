/**
 * @file   ntn_initiation_function.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Mon Sep  2 14:35:11 2013
 *
 * @brief  implementation of initializing ntn and ntrf friction
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
#include "ntn_initiation_function.hh"
#include "ntrf_friction.hh"

// friction regularisations
#include "ntn_fricreg_rubin_ampuero.hh"

// friction laws
#include "ntn_friclaw_linear_slip_weakening.hh"

__BEGIN_SIMTOOLS__

NTNBaseFriction * initializeNTRFFriction(NTRFContact & contact, 
					 ParameterReader & data) {
  AKANTU_DEBUG_IN();
  
  const std::string & friction_law = data.get<std::string>("friction_law");
  const std::string & friction_reg = data.get<std::string>("friction_regularisation");
  
  NTNBaseFriction * friction;
  
  if (friction_law == "coulomb") {
    if (friction_reg == "no_regularisation") {
      friction = new NTRFFriction<NTNFricLawCoulomb, NTNFricRegNoRegularisation>(contact);
    }
    else if (friction_reg == "rubin_ampuero") {
      friction = new NTRFFriction<NTNFricLawCoulomb, NTNFricRegRubinAmpuero>(contact);
      
      friction->setParam("t_star", data.get<Real>("t_star"));
    }
    else {
      AKANTU_DEBUG_ERROR("Do not know the following friction regularisation: " 
			 << friction_reg);
    }

    friction->setParam("mu_s", data.get<Real>("mu_s"));
  }
  else if (friction_law == "linear_slip_weakening") {
    if (friction_reg == "no_regularisation") {
      friction = new NTRFFriction<NTNFricLawLinearSlipWeakening, 
				  NTNFricRegNoRegularisation>(contact);
    }
    else if (friction_reg == "rubin_ampuero") {
      friction = new NTRFFriction<NTNFricLawLinearSlipWeakening, 
				  NTNFricRegRubinAmpuero>(contact);

      friction->setParam("t_star", data.get<Real>("t_star"));
    }
    else {
      AKANTU_DEBUG_ERROR("Do not know the following friction regularisation: " 
			 << friction_reg);
    }

    friction->setParam("mu_s", data.get<Real>("mu_s"));
    friction->setParam("mu_k", data.get<Real>("mu_k"));
    friction->setParam("d_c",  data.get<Real>("d_c"));
  }
  else {
    AKANTU_DEBUG_ERROR("Do not know the following friction law: " 
		       << friction_law);
  }

  AKANTU_DEBUG_OUT();
  return friction;
}

__END_SIMTOOLS__
