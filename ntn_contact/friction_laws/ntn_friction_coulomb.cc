/**
 * @file   ntn_friction_coulomb.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Tue Nov 20 14:23:57 2012
 *
 * @brief  
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
#include "ntn_friction_coulomb.hh"

__BEGIN_SIMTOOLS__

/* -------------------------------------------------------------------------- */
NTNFrictionCoulomb::NTNFrictionCoulomb(NTNContact & contact,
				       const FrictionID & id,
				       const MemoryID & memory_id) : 
  NTNFriction(contact,id,memory_id),
  mu(0,1,0.,id+":mu",0.,"mu") {
  AKANTU_DEBUG_IN();

  NTNFriction::registerSynchronizedArray(this->mu);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFrictionCoulomb::computeFrictionalStrength() {
  AKANTU_DEBUG_IN();

  SolidMechanicsModel & model = this->contact.getModel();
  UInt dim = model.getSpatialDimension();

  UInt nb_ntn_pairs = this->contact.getNbContactNodes();

  // get contact arrays
  const SynchronizedArray<bool> & is_in_contact = this->contact.getIsInContact();
  Real * contact_pressure = this->contact.getContactPressure().storage();

  for (UInt n=0; n<nb_ntn_pairs; ++n) {
    // node pair is NOT in contact
    if (!is_in_contact(n))
      this->frictional_strength(n) = 0.;

    // node pair is in contact
    else {
      // compute frictional strength
      this->frictional_strength(n) = this->mu(n) * Math::norm(dim, &(contact_pressure[n*dim]));
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFrictionCoulomb::registerSynchronizedArray(SynchronizedArrayBase & array) {
  AKANTU_DEBUG_IN();
  
  this->mu.registerDependingArray(array);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFrictionCoulomb::dumpRestart(const std::string & file_name) const {
  AKANTU_DEBUG_IN();
  
  this->mu.dumpRestartFile(file_name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFrictionCoulomb::readRestart(const std::string & file_name) {
  AKANTU_DEBUG_IN();
  
  this->mu.readRestartFile(file_name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFrictionCoulomb::setMu(Real mu) {
  AKANTU_DEBUG_IN();

  NTNFriction::setInternalArray(this->mu, mu);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFrictionCoulomb::setMu(UInt node, Real mu) {
  AKANTU_DEBUG_IN();

  NTNFriction::setInternalArray(this->mu, node, mu);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFrictionCoulomb::printself(std::ostream & stream, int indent) const {
  AKANTU_DEBUG_IN();
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);
  
  stream << space << "NTNFrictionCoulomb [" << std::endl;

  stream << space << this->mu << std::endl;

  stream << space << "]" << std::endl;


  AKANTU_DEBUG_OUT();
}

__END_SIMTOOLS__
