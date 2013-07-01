/**
 * @file   ntrf_friction_mathilde.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Thu May 23 16:34:58 2013
 *
 * @brief  implementation of ntrf friction regularized coulomb
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
#include "ntrf_friction_mathilde.hh"

__BEGIN_SIMTOOLS__

/* -------------------------------------------------------------------------- */
NTRFFrictionMathilde::NTRFFrictionMathilde(NTRFContact & contact,
					   const FrictionID & id,
					   const MemoryID & memory_id) :
  NTRFFriction(contact,id,memory_id),
  t_star(0,1,0.,id+":t_star",0.,"t_star"),
  mu(0,1,0.,id+":mu",0.,"mu"),
  frictional_contact_pressure(0,1,0.,id+":frictional_contact_pressure",0.,
			      "frictionl_contact_pressure"),
  weakening_length(0,1,0.,id+":weakening_lingth",0.,"weakening_length"),
  mu_s(0,1,0.,id+":mu_s",0.,"mu_s"),
  mu_k(0,1,0.,id+":mu_k",0.,"mu_k")
{
  AKANTU_DEBUG_IN();

  NTRFFriction::registerSynchronizedArray(this->t_star);
  NTRFFriction::registerSynchronizedArray(this->weakening_length);
  NTRFFriction::registerSynchronizedArray(this->mu_s);
  NTRFFriction::registerSynchronizedArray(this->mu_k);
  NTRFFriction::registerSynchronizedArray(this->mu);
  NTRFFriction::registerSynchronizedArray(this->frictional_contact_pressure);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFrictionMathilde::computeFrictionalContactPressure() {
  AKANTU_DEBUG_IN();
  
  SolidMechanicsModel & model = this->contact.getModel();
  UInt dim = model.getSpatialDimension();

  UInt nb_contact_nodes = this->contact.getNbContactNodes();

  // get contact arrays
  const SynchronizedArray<bool> & is_in_contact = this->contact.getIsInContact();
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
void NTRFFrictionMathilde::computeFrictionCoefficient() {
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
	// mu = mu_k + (1 - slip / Dc) * (mu_s - mu_k)
	this->mu(n) = this->mu_k(n) 
	            + (1 - (this->slip(n) / this->weakening_length(n))) * (this->mu_s(n) - this->mu_k(n));
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFrictionMathilde::computeFrictionalStrength() {
  this->computeFrictionCoefficient();
  this->computeFrictionalContactPressure();

  SolidMechanicsModel & model = this->contact.getModel();
  UInt dim = model.getSpatialDimension();
  Real delta_t = model.getTimeStep();
  
  UInt nb_contact_nodes = this->contact.getNbContactNodes();

  // get contact arrays
  const SynchronizedArray<bool> & is_in_contact = this->contact.getIsInContact();

  for (UInt n=0; n<nb_contact_nodes; ++n) {
    // node pair is NOT in contact
    if (!is_in_contact(n))
      this->frictional_strength(n) = 0.;

    else {
      // compute frictional contact pressure
      // backward euler method: first order implicit numerical integration method
      // \reg_pres_n+1 = (\reg_pres_n + \delta_t / \t_star * \cur_pres)
      //               / (1 + \delta_t / \t_star)
      Real current_frictional_strength = this->mu(n) * this->frictional_contact_pressure(n);
      Real alpha = delta_t / this->t_star(n);
      this->frictional_strength(n) += alpha * current_frictional_strength;
      this->frictional_strength(n) /= 1 + alpha;
    }
  }
}

/* -------------------------------------------------------------------------- */
void NTRFFrictionMathilde::registerSynchronizedArray(SynchronizedArrayBase & array) {
  AKANTU_DEBUG_IN();

  this->t_star.registerDependingArray(array);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFrictionMathilde::dumpRestart(const std::string & file_name) const {
  AKANTU_DEBUG_IN();

  this->t_star.dumpRestartFile(file_name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFrictionMathilde::readRestart(const std::string & file_name) {
  AKANTU_DEBUG_IN();

  this->t_star.readRestartFile(file_name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFrictionMathilde::setMuS(Real mu) {
  AKANTU_DEBUG_IN();

  NTRFFriction::setInternalArray(this->mu_s, mu);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFrictionMathilde::setMuS(UInt node, Real mu) {
  AKANTU_DEBUG_IN();

  NTRFFriction::setInternalArray(this->mu_s, node, mu);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFrictionMathilde::setMuK(Real mu) {
  AKANTU_DEBUG_IN();

  NTRFFriction::setInternalArray(this->mu_k, mu);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFrictionMathilde::setMuK(UInt node, Real mu) {
  AKANTU_DEBUG_IN();

  NTRFFriction::setInternalArray(this->mu_k, node, mu);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFrictionMathilde::setWeakeningLength(Real length) {
  AKANTU_DEBUG_IN();

  NTRFFriction::setInternalArray(this->weakening_length, length);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFrictionMathilde::setWeakeningLength(UInt node, Real length) {
  AKANTU_DEBUG_IN();

  NTRFFriction::setInternalArray(this->weakening_length, node, length);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFrictionMathilde::setTStar(Real tstar) {
  AKANTU_DEBUG_IN();

  NTRFFriction::setInternalArray(this->t_star, tstar);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFrictionMathilde::setTStar(UInt node, Real tstar) {
  AKANTU_DEBUG_IN();

  NTRFFriction::setInternalArray(this->t_star, node, tstar);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFrictionMathilde::printself(std::ostream & stream, int indent) const {
  AKANTU_DEBUG_IN();
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "NTRFFrictionMathilde [" << std::endl;

  stream << space << this->t_star << std::endl;

  stream << space << "]" << std::endl;


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFrictionMathilde::setToSteadyState() {
  AKANTU_DEBUG_IN();

  this->computeFrictionalContactPressure();
  this->computeFrictionCoefficient();

  SolidMechanicsModel & model = this->contact.getModel();
  UInt dim = model.getSpatialDimension();
  Real delta_t = model.getTimeStep();
  
  UInt nb_contact_nodes = this->contact.getNbContactNodes();

  // get contact arrays
  const SynchronizedArray<bool> & is_in_contact = this->contact.getIsInContact();

  for (UInt n=0; n<nb_contact_nodes; ++n) {
    // node pair is NOT in contact
    if (!is_in_contact(n))
      this->frictional_strength(n) = 0.;

    else {
      this->frictional_strength(n) = this->mu(n) * this->frictional_contact_pressure(n);
    }
  }  
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFrictionMathilde::addDumpFieldToDumper(const std::string & dumper_name,
							  const std::string & field_id) {
  AKANTU_DEBUG_IN();

#ifdef AKANTU_USE_IOHELPER
  //  const SynchronizedArray<UInt> * nodal_filter = &(this->contact.getSlaves());

  if(field_id == "t_star") {
    this->internalAddDumpFieldToDumper(dumper_name,
				       field_id,
				       new DumperIOHelper::NodalField<Real>(this->t_star.getArray()));
  }
  else if(field_id == "mu") {
    this->internalAddDumpFieldToDumper(dumper_name,
				       field_id,
				       new DumperIOHelper::NodalField<Real>(this->mu.getArray()));
  }
  else if (field_id == "frictional_contact_pressure") {
    this->internalAddDumpFieldToDumper(dumper_name,
				       field_id,
				       new DumperIOHelper::NodalField<Real>(this->frictional_contact_pressure.getArray()));
  }
  else if(field_id == "mu_static") {
    this->internalAddDumpFieldToDumper(dumper_name,
				       field_id,
				       new DumperIOHelper::NodalField<Real>(this->mu_s.getArray()));
  }
  else if(field_id == "mu_kinetic") {
    this->internalAddDumpFieldToDumper(dumper_name,
				       field_id,
				       new DumperIOHelper::NodalField<Real>(this->mu_k.getArray()));
  }
  else if(field_id == "weakening_length") {
    this->internalAddDumpFieldToDumper(dumper_name,
				       field_id,
				       new DumperIOHelper::NodalField<Real>(this->weakening_length.getArray()));
  }
  else {
    NTRFFriction::addDumpFieldToDumper(dumper_name, field_id);
  }

#endif

  AKANTU_DEBUG_OUT();
}

__END_SIMTOOLS__
