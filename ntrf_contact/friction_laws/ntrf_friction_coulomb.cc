/**
 * @file   ntrf_friction_coulomb.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Thu Mar 14 14:36:35 2013
 *
 * @brief  implementation of ntrf coulomb friction
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
#include "ntrf_friction_coulomb.hh"

__BEGIN_SIMTOOLS__

/* -------------------------------------------------------------------------- */
NTRFFrictionCoulomb::NTRFFrictionCoulomb(NTRFContact & contact,
					 const FrictionID & id,
					 const MemoryID & memory_id) :
  NTRFFriction(contact,id,memory_id),
  mu(0,1,0.,id+":mu",0.,"mu"),
  frictional_contact_pressure(0,1,0.,id+":frictional_contact_pressure",0.,
			      "frictionl_contact_pressure") {
  AKANTU_DEBUG_IN();

  NTRFFriction::registerSyncronizedArray(this->mu);
  NTRFFriction::registerSyncronizedArray(this->frictional_contact_pressure);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFrictionCoulomb::computeFrictionalContactPressure() {
  AKANTU_DEBUG_IN();
  
  SolidMechanicsModel & model = this->contact.getModel();
  UInt dim = model.getSpatialDimension();

  UInt nb_contact_nodes = this->contact.getNbContactNodes();

  // get contact arrays
  const SyncronizedArray<bool> & is_in_contact = this->contact.getIsInContact();
  Real * contact_pressure = this->contact.getContactPressure().storage();

  for (UInt n=0; n<nb_contact_nodes; ++n) {
    // node pair is NOT in contact
    if (!is_in_contact(n))
      this->frictional_contact_pressure(n) = 0.;

    // node pair is in contact
    else {
      // compute frictional contact pressure
      this->frictional_contact_pressure(n) = Math::norm(dim, &(contact_pressure[n*dim]));
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFrictionCoulomb::computeFrictionalStrength() {
  AKANTU_DEBUG_IN();
  
  this->computeFrictionalContactPressure();

  SolidMechanicsModel & model = this->contact.getModel();
  UInt dim = model.getSpatialDimension();

  UInt nb_contact_nodes = this->contact.getNbContactNodes();

  // get contact arrays
  const SyncronizedArray<bool> & is_in_contact = this->contact.getIsInContact();

  for (UInt n=0; n<nb_contact_nodes; ++n) {
    // node pair is NOT in contact
    if (!is_in_contact(n))
      this->frictional_strength(n) = 0.;

    // node pair is in contact
    else {
      // compute frictional strength
      this->frictional_strength(n) = this->mu(n) * this->frictional_contact_pressure(n);
    }
  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void NTRFFrictionCoulomb::registerSyncronizedArray(SyncronizedArrayBase & array) {
  AKANTU_DEBUG_IN();
  
  this->mu.registerDependingArray(array);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFrictionCoulomb::dumpRestart(const std::string & file_name) const {
  AKANTU_DEBUG_IN();
  
  this->mu.dumpRestartFile(file_name);
  this->frictional_contact_pressure.dumpRestartFile(file_name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFrictionCoulomb::readRestart(const std::string & file_name) {
  AKANTU_DEBUG_IN();
  
  this->mu.readRestartFile(file_name);
  this->frictional_contact_pressure.readRestartFile(file_name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFrictionCoulomb::setMu(Real mu) {
  AKANTU_DEBUG_IN();

  NTRFFriction::setInternalArray(this->mu, mu);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFrictionCoulomb::setMu(UInt node, Real mu) {
  AKANTU_DEBUG_IN();

  NTRFFriction::setInternalArray(this->mu, node, mu);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFrictionCoulomb::printself(std::ostream & stream, int indent) const {
  AKANTU_DEBUG_IN();
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);
  
  stream << space << "NTRFFrictionCoulomb [" << std::endl;

  stream << space << this->mu << std::endl;
  stream << space << this->frictional_contact_pressure << std::endl;

  stream << space << "]" << std::endl;


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFrictionCoulomb::addDumpFieldToDumper(const std::string & dumper_name,
					       const std::string & field_id) {
  AKANTU_DEBUG_IN();
  
#ifdef AKANTU_USE_IOHELPER
  //  const SyncronizedArray<UInt> * nodal_filter = &(this->contact.getSlaves());
  
  if(field_id == "mu") {
    this->internalAddDumpFieldToDumper(dumper_name,
				       field_id,
				       new DumperIOHelper::NodalField<Real>(this->mu.getArray()));
  }
  else if (field_id == "frictional_contact_pressure") {
    this->internalAddDumpFieldToDumper(dumper_name,
				       field_id,
				       new DumperIOHelper::NodalField<Real>(this->frictional_contact_pressure.getArray()));
  }
  else {
    NTRFFriction::addDumpFieldToDumper(dumper_name, field_id);
  }
  
#endif

  AKANTU_DEBUG_OUT();
}

__END_SIMTOOLS__
