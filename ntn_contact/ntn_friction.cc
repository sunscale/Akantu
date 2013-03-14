/**
 * @file   ntn_friction.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Fri Feb 24 15:21:44 2012
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
#include "ntn_friction.hh"

__BEGIN_SIMTOOLS__

/* -------------------------------------------------------------------------- */
NTNFriction::NTNFriction(NTNContact & contact,
			 const FrictionID & id,
			 const MemoryID & memory_id) : 
  Memory(memory_id), id(id),
  contact(contact),
  is_sticking(0,1,true,id+":is_sticking",true,"is_sticking"),
  frictional_strength(0,1,0.,id+":frictional_strength",0.,"frictional_strength"),
  friction_traction(0,contact.getModel().getSpatialDimension(),
		    0.,id+":friction_traction",0.,"friction_traction") {
  //mu(0,1,0.,id+":mu",0.,"mu") {
  AKANTU_DEBUG_IN();

  this->contact.registerSyncronizedArray(this->is_sticking);
  this->contact.registerSyncronizedArray(this->frictional_strength);
  this->contact.registerSyncronizedArray(this->friction_traction);
  //  this->contact.registerSyncronizedArray(this->mu);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFriction::computeFrictionTraction() {
  AKANTU_DEBUG_IN();

  this->computeStickTraction();
  this->computeFrictionalStrength();
  
  SolidMechanicsModel & model = this->contact.getModel();
  UInt dim = model.getSpatialDimension();
  
  UInt nb_ntn_pairs = this->contact.getNbContactNodes();
  
  // get contact arrays
  const SyncronizedArray<bool> & is_in_contact = this->contact.getIsInContact();
  
  // compute friction traction to stop sliding
  Real * contact_pressure = this->contact.getContactPressure().storage();
  Real * friction_trac_p = this->friction_traction.storage();
  this->is_sticking.clear(); // set to not sticking
  
  for (UInt n=0; n<nb_ntn_pairs; ++n) {
    // node pair is NOT in contact
    if (is_in_contact(n)) {
      // check if it is larger than frictional strength
      Real abs_fric = Math::norm(dim, &(friction_trac_p[n*dim]));
      if (abs_fric != 0.) {
	Real alpha = this->frictional_strength(n) / abs_fric;
	
	// larger -> sliding
	if (alpha < 1.) {
	  for (UInt d=0; d<dim; ++d)
	    this->friction_traction(n,d) *= alpha;
	}
	else
	  this->is_sticking(n) = true;
      }
      else {
	// frictional traction is already zero
	this->is_sticking(n) = true;
      }
    }
  }
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFriction::computeStickTraction() {
  AKANTU_DEBUG_IN();

  SolidMechanicsModel & model = this->contact.getModel();
  UInt dim = model.getSpatialDimension();
  Real delta_t = model.getTimeStep();

  UInt nb_ntn_pairs = this->contact.getNbContactNodes();

  // get contact arrays
  const SyncronizedArray<Real> & impedance = this->contact.getImpedance();
  const SyncronizedArray<bool> & is_in_contact = this->contact.getIsInContact();

  // pre-compute the acceleration 
  // (not increment acceleration, because residual is still Kf)
  Array<Real> acceleration(model.getFEM().getMesh().getNbNodes(),dim);
  model.solveLumped(acceleration,
		    model.getMass(),
		    model.getResidual(),
		    model.getBoundary(),
		    model.getF_M2A());

  // compute relative normal fields of velocity and acceleration 
  Array<Real> r_velo(0,dim);
  Array<Real> r_acce(0,dim);
  Array<Real> r_old_acce(0,dim);
  this->contact.computeRelativeTangentialField(model.getVelocity(),     r_velo);
  this->contact.computeRelativeTangentialField(acceleration,            r_acce);
  this->contact.computeRelativeTangentialField(model.getAcceleration(), r_old_acce);

  AKANTU_DEBUG_ASSERT(r_velo.getSize() == nb_ntn_pairs, 
		      "computeRelativeNormalField does not give back arrays " 
		      << "size == nb_ntn_pairs. nb_ntn_pairs = " << nb_ntn_pairs
		      << " | array size = " << r_velo.getSize());

  // compute tangential gap array for all nodes
  Array<Real> gap(nb_ntn_pairs, dim);
  Real * gap_p = gap.storage();
  Real * r_velo_p = r_velo.storage();
  Real * r_acce_p = r_acce.storage();
  Real * r_old_acce_p = r_old_acce.storage();
  for (UInt i=0; i<nb_ntn_pairs*dim; ++i) {
    *gap_p = *r_velo_p + delta_t * *r_acce_p - 0.5 * delta_t * *r_old_acce_p;
    // increment pointers
    gap_p++;
    r_velo_p++;
    r_acce_p++;
    r_old_acce_p++;
  }

  // compute friction traction to stop sliding
  Real * friction_trac_p = this->friction_traction.storage();
  for (UInt n=0; n<nb_ntn_pairs; ++n) {
    
    // node pair is NOT in contact
    if (!is_in_contact(n)) {
      for (UInt d=0; d<dim; ++d)
	this->friction_traction(n,d) = 0.;
    }

    // node pair is in contact
    else {
      // compute friction traction
      for (UInt d=0; d<dim; ++d)
	this->friction_traction(n,d) = impedance(n) * gap(n,d) / 2.;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFriction::applyFrictionTraction() {
  AKANTU_DEBUG_IN();

  SolidMechanicsModel & model = this->contact.getModel();
  Array<Real> & residual = model.getResidual();
  UInt dim = model.getSpatialDimension();

  UInt nb_ntn_pairs = this->contact.getNbContactNodes();
  
  const SyncronizedArray<UInt> & masters = this->contact.getMasters();
  const SyncronizedArray<UInt> & slaves = this->contact.getSlaves();
  const SyncronizedArray<Real> & lumped_boundary = this->contact.getLumpedBoundary();  

  for (UInt n=0; n<nb_ntn_pairs; ++n) {
    UInt master = masters(n);
    UInt slave = slaves(n);
    
    for (UInt d=0; d<dim; ++d) {
      residual(master,d) += lumped_boundary(n,0) * this->friction_traction(n,d);
      residual(slave, d) -= lumped_boundary(n,1) * this->friction_traction(n,d);
    }
  }


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFriction::registerSyncronizedArray(SyncronizedArrayBase & array) {
  AKANTU_DEBUG_IN();
  
  this->frictional_strength.registerDependingArray(array);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFriction::dumpRestart(std::string file_name) const {
  AKANTU_DEBUG_IN();
  
  this->is_sticking.dumpRestartFile(file_name);
  this->frictional_strength.dumpRestartFile(file_name);
  this->friction_traction.dumpRestartFile(file_name);
  //  this->mu.dumpRestartFile(file_name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFriction::readRestart(std::string file_name) {
  AKANTU_DEBUG_IN();
  
  this->is_sticking.readRestartFile(file_name);
  this->frictional_strength.readRestartFile(file_name);
  this->friction_traction.readRestartFile(file_name);
  //  this->mu.readRestartFile(file_name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFriction::setInternalArray(SyncronizedArray<Real> & array, 
				    Real value) {
  AKANTU_DEBUG_IN();

  Real * array_p = array.storage();
  Real nb_array  = array.getSize();
  std::fill_n(array_p, nb_array, value);
  array.setDefaultValue(value);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFriction::setInternalArray(SyncronizedArray<Real> & array, 
				    UInt node, 
				    Real value) {
  AKANTU_DEBUG_IN();

  Int index = this->contact.getNodeIndex(node);
  if (index < 0) {
    AKANTU_DEBUG_WARNING("Node is node a contact node. Therefore, cannot set Mu!!");
  }
  else {
    array(index) = value;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFriction::printself(std::ostream & stream, int indent) const {
  AKANTU_DEBUG_IN();
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);
  
  stream << space << "NTNFriction [" << std::endl;

  stream << space << "]" << std::endl;


  AKANTU_DEBUG_OUT();
}

__END_SIMTOOLS__
