/**
 * @file   ntrf_friction_linear_slip_weakening.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Mon May 13 17:12:09 2013
 *
 * @brief  implementation of linear slip weakening ntrf friction
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
#include "ntrf_friction_linear_slip_weakening.hh"

__BEGIN_SIMTOOLS__

/* -------------------------------------------------------------------------- */
NTRFFrictionLinearSlipWeakening::NTRFFrictionLinearSlipWeakening(NTRFContact & contact,
								 const FrictionID & id,
								 const MemoryID & memory_id) :
  NTRFFrictionRegularizedCoulomb(contact,id,memory_id),
  weakening_length(0,1,0.,id+":weakening_lingth",0.,"weakening_length"),
  mu_s(0,1,0.,id+":mu_s",0.,"mu_s"),
  mu_k(0,1,0.,id+":mu_k",0.,"mu_k") {
  AKANTU_DEBUG_IN();
  
  NTRFFrictionRegularizedCoulomb::registerSyncronizedArray(this->weakening_length);
  NTRFFrictionRegularizedCoulomb::registerSyncronizedArray(this->mu_s);
  NTRFFrictionRegularizedCoulomb::registerSyncronizedArray(this->mu_k);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFrictionLinearSlipWeakening::computeFrictionCoefficient() {
  AKANTU_DEBUG_IN();
  
  SolidMechanicsModel & model = this->contact.getModel();
  UInt dim = model.getSpatialDimension();
  UInt nb_ntn_pairs = this->contact.getNbContactNodes();

  for (UInt n=0; n<nb_ntn_pairs; ++n) {
    if (this->is_sticking(n)) {
      this->mu(n) = this->mu_s(n);
    }
    else {
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
void NTRFFrictionLinearSlipWeakening::computeFrictionalStrength() {
  AKANTU_DEBUG_IN();
  
  computeFrictionCoefficient();
  NTRFFrictionRegularizedCoulomb::computeFrictionalStrength();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFrictionLinearSlipWeakening::registerSyncronizedArray(SyncronizedArrayBase & array) {
  AKANTU_DEBUG_IN();
  
  this->weakening_length.registerDependingArray(array);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFrictionLinearSlipWeakening::dumpRestart(const std::string & file_name) const {
  AKANTU_DEBUG_IN();
  
  this->weakening_length.dumpRestartFile(file_name);
  this->mu_s.dumpRestartFile(file_name);
  this->mu_k.dumpRestartFile(file_name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFrictionLinearSlipWeakening::readRestart(const std::string & file_name) {
  AKANTU_DEBUG_IN();

  this->weakening_length.readRestartFile(file_name);
  this->mu_s.readRestartFile(file_name);
  this->mu_k.readRestartFile(file_name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFrictionLinearSlipWeakening::setMuS(Real mu) {
  AKANTU_DEBUG_IN();

  NTRFFriction::setInternalArray(this->mu_s, mu);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFrictionLinearSlipWeakening::setMuS(UInt node, Real mu) {
  AKANTU_DEBUG_IN();

  NTRFFriction::setInternalArray(this->mu_s, node, mu);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFrictionLinearSlipWeakening::setMuK(Real mu) {
  AKANTU_DEBUG_IN();

  NTRFFriction::setInternalArray(this->mu_k, mu);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFrictionLinearSlipWeakening::setMuK(UInt node, Real mu) {
  AKANTU_DEBUG_IN();

  NTRFFriction::setInternalArray(this->mu_k, node, mu);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFrictionLinearSlipWeakening::setWeakeningLength(Real length) {
  AKANTU_DEBUG_IN();

  NTRFFriction::setInternalArray(this->weakening_length, length);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFrictionLinearSlipWeakening::setWeakeningLength(UInt node, Real length) {
  AKANTU_DEBUG_IN();

  NTRFFriction::setInternalArray(this->weakening_length, node, length);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFrictionLinearSlipWeakening::printself(std::ostream & stream, int indent) const {
  AKANTU_DEBUG_IN();
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);
  
  stream << space << "NTRFFrictionLinearSlipWeakening [" << std::endl;
  NTRFFrictionRegularizedCoulomb::printself(stream, indent+1);

  stream << space << this->weakening_length << std::endl;
  stream << space << this->mu_s << std::endl;
  stream << space << this->mu_k << std::endl;

  stream << space << "]" << std::endl;


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFrictionLinearSlipWeakening::addDumpField(const std::string & field_id) {
  AKANTU_DEBUG_IN();
  
#ifdef AKANTU_USE_IOHELPER
  //  const SyncronizedArray<UInt> * nodal_filter = &(this->contact.getSlaves());
  
  if(field_id == "mu_static") {
    this->contact.addDumpFieldExternal(field_id,
				       new DumperIOHelper::NodalField<Real>(this->mu_s.getArray()));
  }
  else if(field_id == "mu_kinetic") {
    this->contact.addDumpFieldExternal(field_id,
				       new DumperIOHelper::NodalField<Real>(this->mu_k.getArray()));
  }
  else if(field_id == "weakening_length") {
    this->contact.addDumpFieldExternal(field_id,
				       new DumperIOHelper::NodalField<Real>(this->weakening_length.getArray()));
  }
  else {
    NTRFFrictionRegularizedCoulomb::addDumpField(field_id);
  }
  
#endif

  AKANTU_DEBUG_OUT();
}

__END_SIMTOOLS__
