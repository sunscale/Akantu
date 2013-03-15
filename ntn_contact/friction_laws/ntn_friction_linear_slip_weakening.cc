/**
 * @file   ntn_friction_linear_slip_weakening.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Tue Nov 20 15:02:59 2012
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
#include "ntn_friction_linear_slip_weakening.hh"

__BEGIN_SIMTOOLS__

/* -------------------------------------------------------------------------- */
NTNFrictionLinearSlipWeakening::NTNFrictionLinearSlipWeakening(NTNContact & contact,
				       const FrictionID & id,
				       const MemoryID & memory_id) : 
  NTNFrictionCoulomb(contact,id,memory_id),
  weakening_length(0,1,0.,id+":weakening_length",0.,"weakening_length"),
  mu_s(0,1,0.,id+":mu_s",0.,"mu_s"),
  mu_k(0,1,0.,id+":mu_k",0.,"mu_k"),
  slip(0,1,0.,id+":slip",0.,"slip") {
  AKANTU_DEBUG_IN();

  NTNFrictionCoulomb::registerSyncronizedArray(this->weakening_length);
  NTNFrictionCoulomb::registerSyncronizedArray(this->mu_s);
  NTNFrictionCoulomb::registerSyncronizedArray(this->mu_k);
  NTNFrictionCoulomb::registerSyncronizedArray(this->slip);

  // set increment flag of solid mechanics model on
  contact.getModel().setIncrementFlagOn();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFrictionLinearSlipWeakening::computeFrictionCoefficient() {
  AKANTU_DEBUG_IN();
  
  SolidMechanicsModel & model = this->contact.getModel();
  UInt dim = model.getSpatialDimension();
  UInt nb_ntn_pairs = this->contact.getNbContactNodes();

  // compute tangential gap
  Array<Real> rel_tan_incr(0,dim);
  this->contact.computeRelativeTangentialField(model.getIncrement(),
					       rel_tan_incr);
  AKANTU_DEBUG_ASSERT(rel_tan_incr.getSize() == nb_ntn_pairs, 
		      "computeRelativeNormalField does not give back arrays " 
		      << "size == nb_ntn_pairs. nb_ntn_pairs = " << nb_ntn_pairs
		      << " | array size = " << rel_tan_incr.getSize());

  Real * rel_tan_incr_p = rel_tan_incr.storage();
  for (UInt n=0; n<nb_ntn_pairs; ++n) {
    if (this->is_sticking(n)) {
      // reset tangential gap to zero
      this->slip(n) = 0.;
      this->mu(n) = this->mu_s(n);
    }
    else {
      Real abs_incr = Math::norm(dim, &(rel_tan_incr_p[n*dim]));
      this->slip(n) += abs_incr;
      
      if (this->slip(n) >= this->weakening_length(n)) {
	this->mu(n) = this->mu_k(n);
      }
      else {
	// mu = mu_k + (1 - delta / Dc) * (mu_s - mu_k)
	Real delta = this->slip(n) - this->weakening_length(n);
	delta /= this->weakening_length(n);
	this->mu(n) = this->mu_k(n) + (1 - delta) * (this->mu_s(n) - this->mu_k(n));
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFrictionLinearSlipWeakening::computeFrictionalStrength() {
  AKANTU_DEBUG_IN();

  computeFrictionCoefficient();
  NTNFrictionCoulomb::computeFrictionalStrength();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFrictionLinearSlipWeakening::registerSyncronizedArray(SyncronizedArrayBase & array) {
  AKANTU_DEBUG_IN();
  
  this->weakening_length.registerDependingArray(array);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFrictionLinearSlipWeakening::dumpRestart(const std::string & file_name) const {
  AKANTU_DEBUG_IN();
  
  this->weakening_length.dumpRestartFile(file_name);
  this->mu_s.dumpRestartFile(file_name);
  this->mu_k.dumpRestartFile(file_name);
  this->slip.dumpRestartFile(file_name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFrictionLinearSlipWeakening::readRestart(const std::string & file_name) {
  AKANTU_DEBUG_IN();

  this->weakening_length.readRestartFile(file_name);
  this->mu_s.readRestartFile(file_name);
  this->mu_k.readRestartFile(file_name);
  this->slip.readRestartFile(file_name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFrictionLinearSlipWeakening::setMuS(Real mu) {
  AKANTU_DEBUG_IN();

  NTNFriction::setInternalArray(this->mu_s, mu);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFrictionLinearSlipWeakening::setMuS(UInt node, Real mu) {
  AKANTU_DEBUG_IN();

  NTNFriction::setInternalArray(this->mu_s, node, mu);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFrictionLinearSlipWeakening::setMuK(Real mu) {
  AKANTU_DEBUG_IN();

  NTNFriction::setInternalArray(this->mu_k, mu);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFrictionLinearSlipWeakening::setMuK(UInt node, Real mu) {
  AKANTU_DEBUG_IN();

  NTNFriction::setInternalArray(this->mu_k, node, mu);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFrictionLinearSlipWeakening::setWeakeningLength(Real length) {
  AKANTU_DEBUG_IN();

  NTNFriction::setInternalArray(this->weakening_length, length);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFrictionLinearSlipWeakening::setWeakeningLength(UInt node, Real length) {
  AKANTU_DEBUG_IN();

  NTNFriction::setInternalArray(this->weakening_length, node, length);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFrictionLinearSlipWeakening::printself(std::ostream & stream, int indent) const {
  AKANTU_DEBUG_IN();
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);
  
  stream << space << "NTNFrictionLinearSlipWeakening [" << std::endl;
  NTNFrictionCoulomb::printself(stream, indent+1);

  stream << space << this->mu_s << std::endl;
  stream << space << this->mu_k << std::endl;

  stream << space << "]" << std::endl;


  AKANTU_DEBUG_OUT();
}

__END_SIMTOOLS__
