/**
 * @file   ntn_base_contact.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Fri Aug 23 11:05:52 2013
 *
 * @brief  implementation of ntn base contact
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
#include "ntn_base_contact.hh"

__BEGIN_SIMTOOLS__

/* -------------------------------------------------------------------------- */
NTNContactSynchElementFilter::NTNContactSynchElementFilter(NTNBaseContact * contact) :
  contact(contact),
  connectivity(contact->getModel().getMesh().getConnectivities()) {
  AKANTU_DEBUG_IN();
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
bool NTNContactSynchElementFilter::operator()(const Element & e) {
  AKANTU_DEBUG_IN();

  ElementType type = e.type;
  UInt element = e.element;
  GhostType ghost_type = e.ghost_type;
  
  // loop over all nodes of this element
  bool need_element = false;
  UInt nb_nodes = Mesh::getNbNodesPerElement(type);
  for (UInt n=0; n<nb_nodes; ++n) {
    UInt nn = this->connectivity(type, ghost_type)(element,n);
    
    // if one nodes is in this contact, we need this element
    if (this->contact->getNodeIndex(nn) >= 0) {
      need_element = true;
      break;
    }
  }

  AKANTU_DEBUG_OUT();
  return need_element;
}

/* -------------------------------------------------------------------------- */
NTNBaseContact::NTNBaseContact(SolidMechanicsModel & model,
			       const ContactID & id,
			       const MemoryID & memory_id) : 
  Memory(memory_id), Dumpable(), id(id), model(model),
  slaves(0,1,0,id+":slaves",std::numeric_limits<UInt>::quiet_NaN(),"slaves"),
  normals(0,model.getSpatialDimension(),0,id+":normals",
	  std::numeric_limits<Real>::quiet_NaN(),"normals"),
  contact_pressure(0,model.getSpatialDimension(),0,id+":contact_pressure",
		   std::numeric_limits<Real>::quiet_NaN(),"contact_pressure"),
  is_in_contact(0,1,false,id+":is_in_contact",false,"is_in_contact"),
  lumped_boundary_slaves(0,1,0,id+":lumped_boundary_slaves",
			 std::numeric_limits<Real>::quiet_NaN(),"lumped_boundary_slaves"),
  impedance(0,1,0,id+":impedance",std::numeric_limits<Real>::quiet_NaN(),"impedance"),
  contact_surfaces(),
  slave_elements("slave_elements", id, memory_id),
  node_to_elements(),
  synch_registry(NULL)
{
  AKANTU_DEBUG_IN();

  FEM & boundary_fem = this->model.getFEMBoundary();
  for (ghost_type_t::iterator gt = ghost_type_t::begin();  gt != ghost_type_t::end(); ++gt) {
    boundary_fem.initShapeFunctions(*gt);
  }

  Mesh & mesh = this->model.getMesh();
  UInt spatial_dimension = this->model.getSpatialDimension();

  mesh.initByElementTypeArray(this->slave_elements,
			      1,
			      spatial_dimension - 1);

  MeshUtils::buildNode2Elements(mesh,
				this->node_to_elements,
				spatial_dimension - 1);

  this->registerDumper<DumperText>("text_all", id, true);
  this->addDumpFilteredMesh(mesh, 
			    slave_elements, 
			    slaves.getArray(), 
			    spatial_dimension - 1, 
			    _not_ghost,
			    _ek_regular);

  // parallelisation
  this->synch_registry = new SynchronizerRegistry(*this);

  NTNContactSynchElementFilter elem_filter(this);
  this->synchronizer = FilteredSynchronizer::
    createFilteredSynchronizer(this->model.getSynchronizer(),
			       elem_filter);
  this->synch_registry->registerSynchronizer(*(this->synchronizer), _gst_smm_uv);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNBaseContact::findBoundaryElements(const Array<UInt> & interface_nodes, 
					  ByElementTypeArray<UInt> & elements) {
  AKANTU_DEBUG_IN();

  UInt nb_interface_nodes = interface_nodes.getSize();
  const ByElementTypeArray<UInt> & connectivity = this->model.getMesh().getConnectivities();
  
  // add connected boundary elements that have all nodes on this contact
  for (UInt i=0; i<nb_interface_nodes; ++i) {
    UInt node = interface_nodes(i);
    
    CSR<Element>::iterator it   = this->node_to_elements.begin(node);
    CSR<Element>::iterator it_e = this->node_to_elements.end(node);
    for (; it != it_e; ++it) { // loop over all elements connected to node
      ElementType type = it->type;
      UInt element = it->element;
      GhostType ghost_type = it->ghost_type;

      UInt nb_nodes = Mesh::getNbNodesPerElement(type);
      UInt nb_found_nodes = 0;
      for (UInt n=0; n<nb_nodes; ++n) {
	UInt nn = connectivity(type,ghost_type)(element,n);
	if (interface_nodes.find(nn) >= 0)
	  nb_found_nodes++;
	else
	  break;
      }
      
      // this is an element between all contact nodes
      // and is not already in the elements 
      if ((nb_found_nodes == nb_nodes) 
	  && (elements(type,ghost_type).find(element) < 0)) {
	elements(type, ghost_type).push_back(element);
      }
    }
  }
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNBaseContact::addSplitNode(UInt node) {
  AKANTU_DEBUG_IN();
  
  UInt dim = this->model.getSpatialDimension();
  
  // add to node arrays
  this->slaves.push_back(node);
  
  // set contact as false
  this->is_in_contact.push_back(false);
  
  // before initializing
  // set contact pressure, normal, lumped_boundary to Nan
  this->contact_pressure.push_back(std::numeric_limits<Real>::quiet_NaN());
  this->impedance.push_back(std::numeric_limits<Real>::quiet_NaN());
  this->lumped_boundary_slaves.push_back(std::numeric_limits<Real>::quiet_NaN());
  
  Real nan_normal[dim];
  for (UInt d=0; d<dim; ++d)
    nan_normal[d] = std::numeric_limits<Real>::quiet_NaN();
  this->normals.push_back(nan_normal);
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNBaseContact::registerSynchronizedArray(SynchronizedArrayBase & array) {
  AKANTU_DEBUG_IN();
  
  this->slaves.registerDependingArray(array);
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNBaseContact::dumpRestart(const std::string & file_name) const {
  AKANTU_DEBUG_IN();
  
  this->slaves.dumpRestartFile(file_name);
  this->normals.dumpRestartFile(file_name);
  this->is_in_contact.dumpRestartFile(file_name);
  this->contact_pressure.dumpRestartFile(file_name);
  this->lumped_boundary_slaves.dumpRestartFile(file_name);
  this->impedance.dumpRestartFile(file_name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNBaseContact::readRestart(const std::string & file_name) {
  AKANTU_DEBUG_IN();
  
  this->slaves.readRestartFile(file_name);
  this->normals.readRestartFile(file_name);
  this->is_in_contact.readRestartFile(file_name);
  this->contact_pressure.readRestartFile(file_name);
  this->lumped_boundary_slaves.readRestartFile(file_name);
  this->impedance.readRestartFile(file_name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
UInt NTNBaseContact::getNbNodesInContact() const {
  AKANTU_DEBUG_IN();
  
  UInt nb_contact = 0;
  
  UInt nb_nodes = this->getNbContactNodes();
  const Mesh & mesh = this->model.getMesh();

  for (UInt n = 0; n < nb_nodes; ++n) {
    bool is_local_node = mesh.isLocalOrMasterNode(this->slaves(n));
    if (is_local_node && this->is_in_contact(n)) {
      nb_contact++;
    }
  }
  
  StaticCommunicator::getStaticCommunicator().allReduce(&nb_contact, 1, _so_sum);
  
  AKANTU_DEBUG_OUT();
  return nb_contact;
}

/* -------------------------------------------------------------------------- */
void NTNBaseContact::updateInternalData() {
  AKANTU_DEBUG_IN();
  
  updateNormals();
  updateLumpedBoundary();
  updateImpedance();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNBaseContact::updateLumpedBoundary() {
  AKANTU_DEBUG_IN();
  
  this->internalUpdateLumpedBoundary(this->slaves.getArray(),
				     this->slave_elements,
				     this->lumped_boundary_slaves);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNBaseContact::internalUpdateLumpedBoundary(const Array<UInt> & nodes,
						  const ByElementTypeArray<UInt> & elements,
						  SynchronizedArray<Real> & boundary) {
  AKANTU_DEBUG_IN();

  // set all values in lumped_boundary to zero
  boundary.clear();
  
  UInt dim = this->model.getSpatialDimension();
  UInt nb_contact_nodes = getNbContactNodes();
  
  const FEM & boundary_fem = this->model.getFEMBoundary();
  
  const Mesh & mesh = this->model.getMesh();
  
  for (ghost_type_t::iterator gt = ghost_type_t::begin();  gt != ghost_type_t::end(); ++gt) {
    Mesh::type_iterator it = mesh.firstType(dim-1, *gt);
    Mesh::type_iterator last = mesh.lastType(dim-1, *gt);
    for (; it != last; ++it) {
      UInt nb_elements = mesh.getNbElement(*it, *gt);
      UInt nb_nodes_per_element = mesh.getNbNodesPerElement(*it);
      const Array<UInt> & connectivity = mesh.getConnectivity(*it, *gt);

      // get shapes and compute integral
      const Array<Real> & shapes = boundary_fem.getShapes(*it, *gt);
      Array<Real> area(nb_elements,nb_nodes_per_element);
      boundary_fem.integrate(shapes,area,nb_nodes_per_element,*it, *gt);

      if (this->contact_surfaces.size() == 0) {
	AKANTU_DEBUG_WARNING("No surfaces in ntn base contact." 
			     << " You have to define the lumped boundary by yourself.");
      }

      Array<UInt>::const_iterator<UInt> elem_it     = (elements)(*it, *gt).begin();
      Array<UInt>::const_iterator<UInt> elem_it_end = (elements)(*it, *gt).end();
      // loop over contact nodes
      for (; elem_it != elem_it_end; ++elem_it) {
        for (UInt q=0; q<nb_nodes_per_element; ++q) {
          UInt node = connectivity(*elem_it,q);
          UInt node_index = nodes.find(node);
          AKANTU_DEBUG_ASSERT(node_index != -1, 
			      "Could not find node " << node 
			      << " in the array!");
          Real area_to_add = area(*elem_it,q);
          boundary(node_index) += area_to_add;
        }
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNBaseContact::computeContactPressure() {
  AKANTU_DEBUG_IN();

  UInt dim = this->model.getSpatialDimension();
  Real delta_t = this->model.getTimeStep();
  UInt nb_contact_nodes = getNbContactNodes();

  AKANTU_DEBUG_ASSERT(delta_t > 0., 
		      "Cannot compute contact pressure if no time step is set");

  // pre-compute the acceleration 
  // (not increment acceleration, because residual is still Kf)
  Array<Real> acceleration(this->model.getMesh().getNbNodes(),dim);
  this->model.solveLumped(acceleration,
			  this->model.getMass(),
			  this->model.getResidual(),
			  this->model.getBoundary(),
			  this->model.getF_M2A());

  const Array<Real> & residual = this->model.getResidual();

  // compute relative normal fields of displacement, velocity and acceleration 
  Array<Real> r_disp(0,1);
  Array<Real> r_velo(0,1);
  Array<Real> r_acce(0,1);
  Array<Real> r_old_acce(0,1);
  computeNormalGap(r_disp);
  //  computeRelativeNormalField(this->model.getCurrentPosition(), r_disp);
  computeRelativeNormalField(this->model.getVelocity(),        r_velo);
  computeRelativeNormalField(acceleration,                     r_acce);
  computeRelativeNormalField(this->model.getAcceleration(),    r_old_acce);

  AKANTU_DEBUG_ASSERT(r_disp.getSize() == nb_contact_nodes, 
		      "computeRelativeNormalField does not give back arrays " 
		      << "size == nb_contact_nodes. nb_contact_nodes = " << nb_contact_nodes
		      << " | array size = " << r_disp.getSize());

  // compute gap array for all nodes
  Array<Real> gap(nb_contact_nodes, 1);
  Real * gap_p = gap.storage();
  Real * r_disp_p = r_disp.storage();
  Real * r_velo_p = r_velo.storage();
  Real * r_acce_p = r_acce.storage();
  Real * r_old_acce_p = r_old_acce.storage();
  for (UInt i=0; i<nb_contact_nodes; ++i) {
    *gap_p = *r_disp_p + delta_t * *r_velo_p + delta_t * delta_t * *r_acce_p - 0.5 * delta_t * delta_t * *r_old_acce_p;
    // increment pointers
    gap_p++;
    r_disp_p++;
    r_velo_p++;
    r_acce_p++;
    r_old_acce_p++;
  }

  // check if gap is negative -> is in contact
  for (UInt n=0; n<nb_contact_nodes; ++n) {
    if (gap(n) <= 0.) {
      for (UInt d=0; d<dim; ++d) {
	this->contact_pressure(n,d) = this->impedance(n) * gap(n) / (2 * delta_t) 
	                            * this->normals(n,d);
      }
      this->is_in_contact(n) = true;
    }
    else {
      for (UInt d=0; d<dim; ++d)
	this->contact_pressure(n,d) = 0.;
      this->is_in_contact(n) = false;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNBaseContact::applyContactPressure() {
  AKANTU_DEBUG_IN();

  UInt nb_contact_nodes = getNbContactNodes();
  UInt dim = this->model.getSpatialDimension();

  Array<Real> & residual = this->model.getResidual();

  for (UInt n=0; n<nb_contact_nodes; ++n) {
    UInt slave = this->slaves(n);
    
    for (UInt d=0; d<dim; ++d) {
      //residual(master,d) += this->lumped_boundary(n,0) * this->contact_pressure(n,d);
      residual(slave, d) -= this->lumped_boundary_slaves(n) * this->contact_pressure(n,d);
    }
  }
  
  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
Int NTNBaseContact::getNodeIndex(UInt node) const {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
  return this->slaves.find(node);
}

/* -------------------------------------------------------------------------- */
void NTNBaseContact::printself(std::ostream & stream, int indent) const {
  AKANTU_DEBUG_IN();
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);
  
  stream << space << "NTNBaseContact [" << std::endl;
  stream << space << " + id            : " << id << std::endl;
  stream << space << " + slaves        : " << std::endl;
  this->slaves.printself(stream, indent + 2);
  stream << space << " + normals       : " << std::endl;
  this->normals.printself(stream, indent + 2);
  stream << space << " + contact_pressure : " << std::endl;
  this->contact_pressure.printself(stream, indent + 2);

  stream << space << "]" << std::endl;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNBaseContact::syncArrays(SyncChoice sync_choice) {
  AKANTU_DEBUG_IN();
  
  this->slaves.syncElements(sync_choice);
  this->normals.syncElements(sync_choice);
  this->is_in_contact.syncElements(sync_choice);
  this->contact_pressure.syncElements(sync_choice);
  this->lumped_boundary_slaves.syncElements(sync_choice);
  this->impedance.syncElements(sync_choice);
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNBaseContact::addDumpFieldToDumper(const std::string & dumper_name,
					  const std::string & field_id) {
  AKANTU_DEBUG_IN();
  
#ifdef AKANTU_USE_IOHELPER
  const Array<UInt> & nodal_filter = this->slaves.getArray();
  
#define ADD_FIELD(field_id, field, type)				\
  internalAddDumpFieldToDumper(dumper_name,				\
			       field_id,				\
			       new DumperIOHelper::NodalField< type, true, \
							       Array<type>, \
							       Array<UInt> >(field, 0, 0, &nodal_filter))

  if(field_id == "displacement") {
    ADD_FIELD(field_id, this->model.getDisplacement(), Real);
  }
  else if(field_id == "mass") {
    ADD_FIELD(field_id, this->model.getMass(), Real);
  }
  else if(field_id == "velocity") {
    ADD_FIELD(field_id, this->model.getVelocity(), Real);
  }
  else if(field_id == "acceleration") {
    ADD_FIELD(field_id, this->model.getAcceleration(), Real);
  }
  else if(field_id == "force") {
    ADD_FIELD(field_id, this->model.getForce(), Real);
  }
  else if(field_id == "residual") {
    ADD_FIELD(field_id, this->model.getResidual(), Real);
  }
  else if(field_id == "boundary") {
    ADD_FIELD(field_id, this->model.getBoundary(), bool);
  }
  else if(field_id == "increment") {
    ADD_FIELD(field_id, this->model.getIncrement(), Real);
  }
  else if(field_id == "normals") {
    internalAddDumpFieldToDumper(dumper_name,
				 field_id,
				 new DumperIOHelper::NodalField<Real>(this->normals.getArray()));
  }
  else if(field_id == "contact_pressure") {
    internalAddDumpFieldToDumper(dumper_name,
				 field_id,
				 new DumperIOHelper::NodalField<Real>(this->contact_pressure.getArray()));
  }
  else if(field_id == "is_in_contact") {
    internalAddDumpFieldToDumper(dumper_name,
				 field_id,
				 new DumperIOHelper::NodalField<bool>(this->is_in_contact.getArray()));
  }
  else if(field_id == "lumped_boundary_slaves") {
    internalAddDumpFieldToDumper(dumper_name,
				 field_id,
				 new DumperIOHelper::NodalField<Real>(this->lumped_boundary_slaves.getArray()));
  }
  else if(field_id == "impedance") {
    internalAddDumpFieldToDumper(dumper_name,
				 field_id,
				 new DumperIOHelper::NodalField<Real>(this->impedance.getArray()));
  }
  else {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

#undef ADD_FIELD
#endif

  AKANTU_DEBUG_OUT();
}

__END_SIMTOOLS__
