/**
 * @file   ntrf_contact.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Thu Mar 14 14:16:07 2013
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
#include "ntrf_contact.hh"

__BEGIN_SIMTOOLS__

/* -------------------------------------------------------------------------- */
NTRFContact::NTRFContact(SolidMechanicsModel & model,
			 const ContactID & id,
			 const MemoryID & memory_id) :
  Memory(memory_id), Dumpable<DumperParaview>(id), id(id), model(model),
  slaves(0,1,0,id+":slaves",std::numeric_limits<UInt>::quiet_NaN(),"slaves"),
  contact_pressure(0,model.getSpatialDimension(),0,id+":contact_pressure",
		   std::numeric_limits<Real>::quiet_NaN(),"contact_pressure"),
  is_in_contact(0,1,false,id+":is_in_contact",false,"is_in_contact"),
  lumped_boundary(0,1,0,id+":lumped_boundary",
		  std::numeric_limits<Real>::quiet_NaN(),"lumped_boundary"),
  reference_point(0,model.getSpatialDimension()),
  normal(0,model.getSpatialDimension()),
  contact_surfaces(),
  elements("elements", id, memory_id),
  node_to_elements()
{
  AKANTU_DEBUG_IN();

  FEM & boundary_fem = this->model.getFEMBoundary();

  for (ghost_type_t::iterator gt = ghost_type_t::begin();  gt != ghost_type_t::end(); ++gt) {
    boundary_fem.initShapeFunctions(*gt);
  }

  Mesh & mesh = this->model.getMesh();
  UInt spatial_dimension = this->model.getSpatialDimension();

  mesh.initByElementTypeArray(this->elements,
			      1,
			      spatial_dimension - 1);

  MeshUtils::buildNode2Elements(mesh,
				this->node_to_elements,
				spatial_dimension - 1);

  addDumpFilteredMesh(mesh,
		      elements,
		      slaves.getArray(),
		      spatial_dimension - 1,
		      _not_ghost,
		      _ek_regular);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFContact::setReferencePoint(Real x, Real y, Real z) {
  AKANTU_DEBUG_IN();

  Real coord[3];
  coord[0] = x;
  coord[1] = y;
  coord[2] = z;

  this->reference_point.resize(1);

  UInt dim = this->model.getSpatialDimension();
  for (UInt d=0; d<dim; ++d)
    this->reference_point(0,d) = coord[d];

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFContact::setNormal(Real x, Real y, Real z) {
  AKANTU_DEBUG_IN();

  UInt dim = this->model.getSpatialDimension();

  Real coord[3];
  coord[0] = x;
  coord[1] = y;
  coord[2] = z;

  if (dim == 2)
    Math::normalize2(coord);
  else if (dim == 3)
    Math::normalize3(coord);
  else
    AKANTU_DEBUG_ERROR("setNormal in NTRFContact not implemented for dimension " << dim);

  this->normal.resize(1);

  for (UInt d=0; d<dim; ++d)
    this->normal(0,d) = coord[d];


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFContact::addSurface(const Surface & surf) {
  AKANTU_DEBUG_IN();

  UInt dim = this->model.getSpatialDimension();

  const Mesh & mesh_ref = this->model.getFEM().getMesh();

  const SubBoundary & boundary = mesh_ref.getSubBoundary(surf);
  this->contact_surfaces.insert(&boundary);

  // find slave nodes
  for(SubBoundary::nodes_const_iterator nodes_it(boundary.nodes_begin());
      nodes_it!= boundary.nodes_end();
      ++nodes_it) {
    this->addNode(*nodes_it);
  }

  // synchronize with depending nodes
  updateInternalData();
  syncArrays(_added);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFContact::addNode(UInt node) {
  AKANTU_DEBUG_IN();

  UInt dim = this->model.getSpatialDimension();

  // add to node arrays
  this->slaves.push_back(node);

  // set contact as false
  this->is_in_contact.push_back(false);

  // before initializing
  // set contact pressure, normal, lumped_boundary to Nan
  this->contact_pressure.push_back(std::numeric_limits<Real>::quiet_NaN());
  this->lumped_boundary.push_back(std::numeric_limits<Real>::quiet_NaN());

  // add connected boundary elements that have all nodes on this contact
  const ByElementTypeArray<UInt> & connectivity = this->model.getMesh().getConnectivities();
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
      if (this->slaves.find(nn) >= 0)
	nb_found_nodes++;
      else
	break;
    }

    // this is an element between all contact nodes
    if (nb_found_nodes == nb_nodes) {
      this->elements(type, ghost_type).push_back(element);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFContact::addNodes(Array<UInt> & nodes) {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = nodes.getSize();
  UInt nb_compo = nodes.getNbComponent();
  for (UInt n=0; n<nb_nodes; ++n) {
    for (UInt c=0; c<nb_compo; ++c) {
      this->addNode(nodes(n,c));
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFContact::registerSyncronizedArray(SyncronizedArrayBase & array) {
  AKANTU_DEBUG_IN();

  this->slaves.registerDependingArray(array);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFContact::dumpRestart(const std::string & file_name) const {
  AKANTU_DEBUG_IN();

  this->slaves.dumpRestartFile(file_name);
  this->is_in_contact.dumpRestartFile(file_name);
  this->contact_pressure.dumpRestartFile(file_name);
  this->lumped_boundary.dumpRestartFile(file_name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFContact::readRestart(const std::string & file_name) {
  AKANTU_DEBUG_IN();

  this->slaves.readRestartFile(file_name);
  this->is_in_contact.readRestartFile(file_name);
  this->contact_pressure.readRestartFile(file_name);
  this->lumped_boundary.readRestartFile(file_name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFContact::updateInternalData() {
  AKANTU_DEBUG_IN();

  updateLumpedBoundary();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFContact::updateLumpedBoundary() {
  AKANTU_DEBUG_IN();

  // set all values in lumped_boundary to zero
  this->lumped_boundary.clear();

  UInt dim = this->model.getSpatialDimension();
  UInt nb_contact_nodes = getNbContactNodes();

  const FEM & boundary_fem = this->model.getFEMBoundary();

  const Mesh & mesh = this->model.getFEM().getMesh();

  // get elements connected to each node
  const CSR<Element> & node_to_element = this->node_to_elements;

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
      debug::setDebugLevel(dblDump);
      area.printself(std::cout);
      debug::setDebugLevel(dblWarning);
      // get surface id information
  //    const Array<UInt> & surface_id = mesh.getSurfaceID(*it);
  //    std::set<UInt>::iterator pos;
  //    std::set<UInt>::iterator end = this->contact_surfaces.end();

      if (this->contact_surfaces.size() == 0)
        std::cerr << "No surfaces in ntrf contact. You have to define the lumped boundary by yourself." << std::endl;

      // loop over contact nodes
      for (UInt i=0; i<nb_contact_nodes; ++i) {
        UInt node = this->slaves(i);

        CSR<Element>::const_iterator elem = node_to_element.begin(node);
        // loop over all elements connected to this node
        for (; elem != node_to_element.end(node); ++elem) {
          Element e = *elem;
          if(e.ghost_type != *gt) {
            continue;
          }
            // if element is not at interface continue
          //	pos = this->contact_surfaces.find(surface_id(e));
          //	if (pos == end)
          //	  continue;

          // loop over all points of this element
          for (UInt q=0; q<nb_nodes_per_element; ++q) {
            if (connectivity(e.element,q) == node) {
              Real area_to_add = area(e.element,q);
              this->lumped_boundary(i) += area_to_add;
            }
          }
        }
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFContact::computeContactPressure() {
  AKANTU_DEBUG_IN();

  UInt dim = this->model.getSpatialDimension();
  Real delta_t = this->model.getTimeStep();
  UInt nb_contact_nodes = getNbContactNodes();

  AKANTU_DEBUG_ASSERT(delta_t > 0.,
		      "Cannot compute contact pressure if no time step is set");

  // pre-compute the acceleration
  // (not increment acceleration, because residual is still Kf)
  Array<Real> acceleration(this->model.getFEM().getMesh().getNbNodes(),dim);
  this->model.solveLumped(acceleration,
			  this->model.getMass(),
			  this->model.getResidual(),
			  this->model.getBoundary(),
			  this->model.getF_M2A());

  const Array<Real> & residual = this->model.getResidual();
  const Array<Real> & mass = this->model.getMass();

  // compute relative normal fields of displacement, velocity and acceleration
  Array<Real> r_disp(0,1);
  Array<Real> r_velo(0,1);
  Array<Real> r_acce(0,1);
  Array<Real> r_old_acce(0,1);
  computeNormalGap(r_disp);
  computeRelativeNormalField(this->model.getVelocity(),        r_velo);
  computeRelativeNormalField(acceleration,                     r_acce);
  computeRelativeNormalField(this->model.getAcceleration(),    r_old_acce);

  AKANTU_DEBUG_ASSERT(r_disp.getSize() == nb_contact_nodes,
		      "computeNormalGap does not give back arrays "
		      << "size == nb_contact_nodes. nb_contact_nodes = "
		      << nb_contact_nodes << " | array size = "
		      << r_disp.getSize());
  AKANTU_DEBUG_ASSERT(r_velo.getSize() == nb_contact_nodes,
		      "computeRelativeNormalField does not give back arrays "
		      << "size == nb_contact_nodes. nb_contact_nodes = "
		      << nb_contact_nodes << " | array size = "
		      << r_velo.getSize());

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
    if (gap(n) < 0.) {
      UInt node = this->slaves(n);
      for (UInt d=0; d<dim; ++d) {
	this->contact_pressure(n,d) = mass(node) / delta_t / delta_t
	  / this->lumped_boundary(n) * gap(n) * this->normal(0,d);
	/*
	// from the ntn_contact:
	this->contact_pressure(n,d) = this->impedance(n) * gap(n) / (2 * delta_t)
	* this->normal(0,d);
	*/
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
void NTRFContact::applyContactPressure() {
  AKANTU_DEBUG_IN();

  UInt nb_contact_nodes = getNbContactNodes();
  UInt dim = this->model.getSpatialDimension();

  Array<Real> & residual = this->model.getResidual();

  for (UInt n=0; n<nb_contact_nodes; ++n) {
    UInt node = this->slaves(n);
    for (UInt d=0; d<dim; ++d) {
      residual(node, d) -= this->lumped_boundary(n) * this->contact_pressure(n,d);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFContact::computeRelativeTangentialField(const Array<Real> & field,
						Array<Real> & rel_tang_field) const {
  AKANTU_DEBUG_IN();

  // resize arrays to zero
  rel_tang_field.resize(0);

  UInt dim = this->model.getSpatialDimension();
  Real * field_p = field.storage();
  Real * normals_p = this->normal.storage();

  UInt nb_contact_nodes = this->slaves.getSize();
  for (UInt n=0; n<nb_contact_nodes; ++n) {
    // nodes
    UInt node = this->slaves(n);

    // compute relative field to master
    Real rel_array[dim];
    for (UInt d=0; d<dim; ++d) {
      rel_array[d] = field_p[node*dim + d];
    }

    // compute dot product with normal of master
    Real dot_prod = Math::vectorDot(normals_p, rel_array, dim);

    // compute the tangential projection of the relative field to the master
    for (UInt d=0; d<dim; ++d)
      rel_array[d] -= dot_prod * normals_p[d];

    rel_tang_field.push_back(rel_array);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFContact::computeNormalGap(Array<Real> & gap) const {
  AKANTU_DEBUG_IN();

  gap.resize(0);

  UInt dim = this->model.getSpatialDimension();
  Real * cur_pos = this->model.getCurrentPosition().storage();
  Real * normals_p = this->normal.storage();

  UInt nb_contact_nodes = this->getNbContactNodes();
  for (UInt n=0; n<nb_contact_nodes; ++n) {
    // nodes
    UInt node = this->slaves(n);

    // compute relative field to master
    Real rel_array[dim];
    for (UInt d=0; d<dim; ++d) {
      rel_array[d] = cur_pos[node*dim + d] - this->reference_point(0,d);
    }

    // compute dot product with normal of master
    Real dot_prod = Math::vectorDot(normals_p, rel_array, dim);

    gap.push_back(dot_prod);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFContact::computeRelativeNormalField(const Array<Real> & field,
					    Array<Real> & rel_normal_field) const {
  AKANTU_DEBUG_IN();

  // resize arrays to zero
  rel_normal_field.resize(0);

  UInt dim = this->model.getSpatialDimension();
  Real * field_p = field.storage();
  Real * normals_p = this->normal.storage();

  UInt nb_contact_nodes = this->getNbContactNodes();
  for (UInt n=0; n<nb_contact_nodes; ++n) {
    // nodes
    UInt node = this->slaves(n);

    // compute dot product with normal of master
    Real dot_prod = Math::vectorDot(normals_p, &(field_p[node*dim]), dim);

    rel_normal_field.push_back(dot_prod);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Int NTRFContact::getNodeIndex(UInt node) const {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
  return this->slaves.find(node);
}

/* -------------------------------------------------------------------------- */
void NTRFContact::printself(std::ostream & stream, int indent) const {
  AKANTU_DEBUG_IN();
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "NTRFContact [" << std::endl;
  stream << space << " + id            : " << id << std::endl;
  stream << space << " + slaves        : " << std::endl;
  this->slaves.printself(stream, indent + 2);
  stream << space << " + contact_pressure : " << std::endl;
  this->contact_pressure.printself(stream, indent + 2);

  stream << space << "]" << std::endl;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFContact::syncArrays(SyncChoice sync_choice) {
  AKANTU_DEBUG_IN();

  this->slaves.syncElements(sync_choice);
  this->is_in_contact.syncElements(sync_choice);
  this->contact_pressure.syncElements(sync_choice);
  this->lumped_boundary.syncElements(sync_choice);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFContact::addDumpField(const std::string & field_id) {
  AKANTU_DEBUG_IN();

#ifdef AKANTU_USE_IOHELPER
  const Array<UInt> & nodal_filter = this->slaves.getArray();

#define ADD_FIELD(field_id, field, type)				\
  addDumpFieldToDumper(field_id,					\
		       new DumperIOHelper::NodalField< type, true,	\
						       Array<type>,	\
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
  else if(field_id == "contact_pressure") {
    addDumpFieldToDumper(field_id,
			 new DumperIOHelper::NodalField<Real>(this->contact_pressure.getArray()));
  }
  else if(field_id == "is_in_contact") {
    addDumpFieldToDumper(field_id,
			 new DumperIOHelper::NodalField<bool>(this->is_in_contact.getArray()));
  }
  else if(field_id == "lumped_boundary") {
    addDumpFieldToDumper(field_id,
			 new DumperIOHelper::NodalField<Real>(this->lumped_boundary.getArray()));
  }
  else {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

#undef ADD_FIELD
#endif

  AKANTU_DEBUG_OUT();
}

__END_SIMTOOLS__
