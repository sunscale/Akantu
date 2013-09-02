/**
 * @file   ntn_friction_tmpl.hh
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

__BEGIN_SIMTOOLS__

/* -------------------------------------------------------------------------- */
NTNFriction::NTNFriction(NTNContact & contact,
			 const FrictionID & id,
			 const MemoryID & memory_id) : 
  NTNBaseFriction(&contact, id, memory_id) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFriction::applyFrictionTraction() {
  AKANTU_DEBUG_IN();

  NTNContact * ntn_contact = dynamic_cast<NTNContact * >(this->contact);
  SolidMechanicsModel & model = ntn_contact->getModel();
  Array<Real> & residual = model.getResidual();
  UInt dim = model.getSpatialDimension();
  
  const SynchronizedArray<UInt> & masters = ntn_contact->getMasters();
  const SynchronizedArray<UInt> & slaves  = ntn_contact->getSlaves();
  const SynchronizedArray<Real> & l_boundary_slaves  = ntn_contact->getLumpedBoundarySlaves();
  const SynchronizedArray<Real> & l_boundary_masters = ntn_contact->getLumpedBoundaryMasters();

  UInt nb_contact_nodes = ntn_contact->getNbContactNodes();
  for (UInt n=0; n<nb_contact_nodes; ++n) {
    UInt master = masters(n);
    UInt slave  = slaves(n);
    
    for (UInt d=0; d<dim; ++d) {
      residual(master,d) += l_boundary_masters(n) * this->friction_traction(n,d);
      residual(slave, d) -= l_boundary_slaves(n)  * this->friction_traction(n,d);
    }
  }


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFriction::printself(std::ostream & stream, int indent) const {
  AKANTU_DEBUG_IN();
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);
  
  stream << space << "NTNFriction [" << std::endl;
  NTNBaseFriction::printself(stream, indent);
  stream << space << "]" << std::endl;


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/*
void NTNFriction::addDumpFieldToDumper(const std::string & dumper_name,
				       const std::string & field_id) {
  AKANTU_DEBUG_IN();
  
#ifdef AKANTU_USE_IOHELPER
  //  const SynchronizedArray<UInt> * nodal_filter = &(this->contact->getSlaves());
  
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
  else if(field_id == "slip_speed") {
    this->internalAddDumpFieldToDumper(dumper_name,
			       field_id,
			       new DumperIOHelper::NodalField<Real>(this->slip_speed.getArray()));
  }
  else {
    this->contact->addDumpFieldToDumper(dumper_name, field_id);
  }
  
#endif

  AKANTU_DEBUG_OUT();
}
*/

__END_SIMTOOLS__
