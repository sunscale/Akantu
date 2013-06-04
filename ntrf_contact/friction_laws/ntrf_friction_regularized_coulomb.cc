/**
 * @file   ntrf_friction_regularized_coulomb.cc
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
#include "ntrf_friction_regularized_coulomb.hh"

__BEGIN_SIMTOOLS__

/* -------------------------------------------------------------------------- */
NTRFFrictionRegularizedCoulomb::NTRFFrictionRegularizedCoulomb(NTRFContact & contact,
							       const FrictionID & id,
							       const MemoryID & memory_id) :
  NTRFFrictionCoulomb(contact,id,memory_id),
  regularization_on(false),
  t_star(0,1,0.,id+":t_star",0.,"t_star") {
  AKANTU_DEBUG_IN();

  NTRFFrictionCoulomb::registerSyncronizedArray(this->t_star);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFrictionRegularizedCoulomb::computeFrictionalContactPressure() {
  AKANTU_DEBUG_IN();

  if (!this->regularization_on)
    NTRFFrictionCoulomb::computeFrictionalContactPressure();
  else {
    SolidMechanicsModel & model = this->contact.getModel();
    UInt dim = model.getSpatialDimension();
    Real delta_t = model.getTimeStep();

    UInt nb_contact_nodes = this->contact.getNbContactNodes();

    // get contact arrays
    const SyncronizedArray<bool> & is_in_contact = this->contact.getIsInContact();
    Real * contact_pressure = this->contact.getContactPressure().storage();

    for (UInt n=0; n<nb_contact_nodes; ++n) {
      // node pair is NOT in contact
      if (!is_in_contact(n))
	this->frictional_contact_pressure(n) = 0.;

      // if t_star is too small compute like Coulomb friction (without regularization)
      else if (Math::are_float_equal(this->t_star(n), 0.)) {
	this->frictional_contact_pressure(n) = Math::norm(dim, &(contact_pressure[n*dim]));
      }

      else {
	// compute frictional contact pressure
	// backward euler method: first order implicit numerical integration method
	// \reg_pres_n+1 = (\reg_pres_n + \delta_t / \t_star * \cur_pres)
	//               / (1 + \delta_t / \t_star)
	Real current_contact_pressure = Math::norm(dim, &(contact_pressure[n*dim]));
	Real alpha = delta_t / this->t_star(n);
	this->frictional_contact_pressure(n) += alpha * current_contact_pressure;
	this->frictional_contact_pressure(n) /= 1 + alpha;
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFrictionRegularizedCoulomb::registerSyncronizedArray(SyncronizedArrayBase & array) {
  AKANTU_DEBUG_IN();

  this->t_star.registerDependingArray(array);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFrictionRegularizedCoulomb::dumpRestart(const std::string & file_name) const {
  AKANTU_DEBUG_IN();

  this->t_star.dumpRestartFile(file_name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFrictionRegularizedCoulomb::readRestart(const std::string & file_name) {
  AKANTU_DEBUG_IN();

  this->t_star.readRestartFile(file_name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFrictionRegularizedCoulomb::setTStar(Real tstar) {
  AKANTU_DEBUG_IN();

  NTRFFriction::setInternalArray(this->t_star, tstar);
  this->regularization_on = true;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFrictionRegularizedCoulomb::setTStar(UInt node, Real tstar) {
  AKANTU_DEBUG_IN();

  NTRFFriction::setInternalArray(this->t_star, node, tstar);
  this->regularization_on = true;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFrictionRegularizedCoulomb::printself(std::ostream & stream, int indent) const {
  AKANTU_DEBUG_IN();
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "NTRFFrictionRegularizedCoulomb [" << std::endl;

  stream << space << "Regularization On: " << this->regularization_on << std::endl;
  stream << space << this->t_star << std::endl;

  stream << space << "]" << std::endl;


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFrictionRegularizedCoulomb::setToSteadyState() {
  AKANTU_DEBUG_IN();
  NTRFFrictionCoulomb::computeFrictionalContactPressure();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFrictionRegularizedCoulomb::addDumpField(const std::string & field_id) {
  AKANTU_DEBUG_IN();

#ifdef AKANTU_USE_IOHELPER
  //  const SyncronizedArray<UInt> * nodal_filter = &(this->contact.getSlaves());

  if(field_id == "t_star") {
    this->contact.addDumpFieldExternal(field_id,
				       new DumperIOHelper::NodalField<Real>(this->t_star.getArray()));
  }
  else {
    NTRFFrictionCoulomb::addDumpField(field_id);
  }

#endif

  AKANTU_DEBUG_OUT();
}

__END_SIMTOOLS__
