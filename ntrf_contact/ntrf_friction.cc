/**
 * @file   ntrf_friction.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Mon Nov  5 10:23:47 2012
 *
 * @brief  implementation of node to rigid flat interface friction
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
#include "ntrf_friction.hh"

__BEGIN_SIMTOOLS__

/* -------------------------------------------------------------------------- */
NTRFFriction::NTRFFriction(NTRFContact & contact,
			   const FrictionID & id,
			   const MemoryID & memory_id) : 
  Memory(memory_id), id(id),
  Dumpable(),
  contact(contact),
  is_sticking(0,1,true,id+":is_sticking",true,"is_sticking"),
  frictional_strength(0,1,0.,id+":frictional_strength",0.,"frictional_strength"),
  friction_traction(0,contact.getModel().getSpatialDimension(),
		    0.,id+":friction_traction",0.,"friction_traction"),
  slip(0,1,0.,id+":slip",0.,"slip") {
  AKANTU_DEBUG_IN();

  this->contact.registerSynchronizedArray(this->is_sticking);
  this->contact.registerSynchronizedArray(this->frictional_strength);
  this->contact.registerSynchronizedArray(this->friction_traction);
  this->contact.registerSynchronizedArray(this->slip);

  contact.getModel().setIncrementFlagOn();

  this->registerExternalDumper(&(contact.getDumper()),
			       contact.getDefaultDumperName(),
			       true);
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFriction::updateSlip() {
  AKANTU_DEBUG_IN();
  
  SolidMechanicsModel & model = this->contact.getModel();
  UInt dim = model.getSpatialDimension();
  UInt nb_nodes = this->contact.getNbContactNodes();
  
  Array<Real> rel_tan_incr(0,dim);
  this->contact.computeRelativeTangentialField(model.getIncrement(),
					       rel_tan_incr);

  Real * rel_tan_incr_p = rel_tan_incr.storage();
  for (UInt n=0; n<nb_nodes; ++n) {
    if (this->is_sticking(n)) {
      this->slip(n) = 0.;
    }
    else {
      Real abs_incr = Math::norm(dim, &(rel_tan_incr_p[n*dim]));
      this->slip(n) += abs_incr;
    }
  }
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFriction::computeFrictionTraction() {
  AKANTU_DEBUG_IN();

  this->computeStickTraction();
  this->computeFrictionalStrength();

  SolidMechanicsModel & model = this->contact.getModel();
  UInt dim = model.getSpatialDimension();

  UInt nb_imp_nodes = this->contact.getNbContactNodes();

  // get contact arrays
  const SynchronizedArray<bool> & is_in_contact = this->contact.getIsInContact();

  // compute friction traction to stop sliding
  Real * contact_pressure = this->contact.getContactPressure().storage();
  Real * friction_trac_p = this->friction_traction.storage();
  this->is_sticking.clear(); // set to not sticking

  for (UInt n=0; n<nb_imp_nodes; ++n) {
    // node pair is in contact
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
void NTRFFriction::computeStickTraction() {
  AKANTU_DEBUG_IN();

  SolidMechanicsModel & model = this->contact.getModel();
  UInt dim = model.getSpatialDimension();
  Real delta_t = model.getTimeStep();

  UInt nb_imp_nodes = this->contact.getNbContactNodes();

  // get model arrays
  const Array<Real> & mass = model.getMass();

  // get contact arrays
  const SynchronizedArray<UInt> & nodes = this->contact.getSlaves();
  const SynchronizedArray<Real> & lumped_boundary = this->contact.getLumpedBoundarySlaves();
  const SynchronizedArray<bool> & is_in_contact = this->contact.getIsInContact();

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

  AKANTU_DEBUG_ASSERT(r_velo.getSize() == nb_imp_nodes, 
		      "computeRelativeNormalField does not give back arrays " 
		      << "size == nb_imp_nodes. nb_imp_nodes = " << nb_imp_nodes
		      << " | array size = " << r_velo.getSize());

  // compute tangential gap array for all nodes
  Array<Real> gap_dot(nb_imp_nodes, dim);
  Real * gap_dot_p = gap_dot.storage();
  Real * r_velo_p = r_velo.storage();
  Real * r_acce_p = r_acce.storage();
  Real * r_old_acce_p = r_old_acce.storage();
  for (UInt i=0; i<nb_imp_nodes*dim; ++i) {
    *gap_dot_p = *r_velo_p + delta_t * *r_acce_p - 0.5 * delta_t * *r_old_acce_p;
    // increment pointers
    gap_dot_p++;
    r_velo_p++;
    r_acce_p++;
    r_old_acce_p++;
  }

  // compute friction traction to stop sliding
  Real * friction_trac_p = this->friction_traction.storage();
  for (UInt n=0; n<nb_imp_nodes; ++n) {
    
    // node pair is NOT in contact
    if (!is_in_contact(n)) {
      for (UInt d=0; d<dim; ++d)
	this->friction_traction(n,d) = 0.;
    }

    // node pair is in contact
    else {
      UInt node = nodes(n);
      // compute friction traction
      for (UInt d=0; d<dim; ++d)
	this->friction_traction(n,d) = mass(node) * gap_dot(n,d) / delta_t / lumped_boundary(n);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFriction::applyFrictionTraction() {
  AKANTU_DEBUG_IN();

  SolidMechanicsModel & model = this->contact.getModel();
  Array<Real> & residual = model.getResidual();
  UInt dim = model.getSpatialDimension();

  UInt nb_imp_nodes = this->contact.getNbContactNodes();
  
  const SynchronizedArray<UInt> & nodes = this->contact.getSlaves();
  const SynchronizedArray<Real> & lumped_boundary = this->contact.getLumpedBoundarySlaves();

  for (UInt n=0; n<nb_imp_nodes; ++n) {
    UInt node = nodes(n);
    for (UInt d=0; d<dim; ++d) {
      residual(node, d) -= lumped_boundary(n) * this->friction_traction(n,d);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFriction::registerSynchronizedArray(SynchronizedArrayBase & array) {
  AKANTU_DEBUG_IN();
  
  this->frictional_strength.registerDependingArray(array);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFriction::dumpRestart(const std::string & file_name) const {
  AKANTU_DEBUG_IN();
  
  this->is_sticking.dumpRestartFile(file_name);
  this->frictional_strength.dumpRestartFile(file_name);
  this->friction_traction.dumpRestartFile(file_name);
  this->slip.dumpRestartFile(file_name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFriction::readRestart(const std::string & file_name) {
  AKANTU_DEBUG_IN();
  
  this->is_sticking.readRestartFile(file_name);
  this->frictional_strength.readRestartFile(file_name);
  this->friction_traction.readRestartFile(file_name);
  this->slip.readRestartFile(file_name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFriction::setInternalArray(SynchronizedArray<Real> & array, 
				    Real value) {
  AKANTU_DEBUG_IN();

  Real * array_p = array.storage();
  Real nb_array  = array.getSize();
  std::fill_n(array_p, nb_array, value);
  array.setDefaultValue(value);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFriction::setInternalArray(SynchronizedArray<Real> & array, 
				    UInt node, 
				    Real value) {
  AKANTU_DEBUG_IN();

  Int index = this->contact.getNodeIndex(node);
  if (index < 0) {
    AKANTU_DEBUG_WARNING("Node is not a contact node. Therefore, cannot set Mu!!");
  }
  else {
    array(index) = value;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
UInt NTRFFriction::getNbStickingNodes() const {
  AKANTU_DEBUG_IN();
  
  UInt nb_stick = 0;

  UInt nb_nodes = this->contact.getNbContactNodes();
  const SynchronizedArray<UInt> & nodes = this->contact.getSlaves();
  const SynchronizedArray<bool> & is_in_contact = this->contact.getIsInContact();

  const Mesh & mesh = this->contact.getModel().getMesh();

  for (UInt n = 0; n < nb_nodes; ++n) {
    bool is_local_node = mesh.isLocalOrMasterNode(nodes(n));
    if (is_local_node && is_in_contact(n) && this->is_sticking(n)) {
      nb_stick++;
    }
  }

  StaticCommunicator::getStaticCommunicator().allReduce(&nb_stick, 1, _so_sum);

  AKANTU_DEBUG_OUT();
  return nb_stick;
}

/* -------------------------------------------------------------------------- */
void NTRFFriction::printself(std::ostream & stream, int indent) const {
  AKANTU_DEBUG_IN();
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);
  
  stream << space << "NTRFFriction [" << std::endl;

  stream << space << "]" << std::endl;


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFFriction::addDumpFieldToDumper(const std::string & dumper_name,
					const std::string & field_id) {
  AKANTU_DEBUG_IN();
  
#ifdef AKANTU_USE_IOHELPER
  //  const SynchronizedArray<UInt> * nodal_filter = &(this->contact.getSlaves());
  
  if(field_id == "is_sticking") {
    this->internalAddDumpFieldToDumper(dumper_name,
			       field_id,
			       new DumperIOHelper::NodalField<bool>(this->is_sticking.getArray()));
  }
  else if(field_id == "frictional_strength") {
    this->internalAddDumpFieldToDumper(dumper_name,
			       field_id,
			       new DumperIOHelper::NodalField<Real>(this->frictional_strength.getArray()));
  }
  else if(field_id == "friction_traction") {
    this->internalAddDumpFieldToDumper(dumper_name,
			       field_id,
			       new DumperIOHelper::NodalField<Real>(this->friction_traction.getArray()));
  }
  else if(field_id == "slip") {
    this->internalAddDumpFieldToDumper(dumper_name,
			       field_id,
			       new DumperIOHelper::NodalField<Real>(this->slip.getArray()));
  }
  else {
    this->contact.addDumpFieldToDumper(dumper_name, field_id);
  }
  
#endif

  AKANTU_DEBUG_OUT();
}

__END_SIMTOOLS__
