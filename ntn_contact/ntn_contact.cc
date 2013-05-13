/**
 * @file   ntn_contact.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Thu Mar 14 11:52:00 2013
 *
 * @brief  implementation of ntn_contact
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
#include "ntn_contact.hh"

__BEGIN_SIMTOOLS__

/* -------------------------------------------------------------------------- */
NTNContact::NTNContact(SolidMechanicsModel & model,
                       const ContactID & id,
                       const MemoryID & memory_id) : 
  Memory(memory_id), id(id), model(model),
  slaves(0,1,0,id+":slaves",std::numeric_limits<UInt>::quiet_NaN(),"slaves"),
  masters(0,1,0,id+":masters",std::numeric_limits<UInt>::quiet_NaN(),"masters"),
  normals(0,model.getSpatialDimension(),0,id+":normals",
  std::numeric_limits<Real>::quiet_NaN(),"normals"),
  contact_pressure(0,model.getSpatialDimension(),0,id+":contact_pressure",
  std::numeric_limits<Real>::quiet_NaN(),"contact_pressure"),
  is_in_contact(0,1,false,id+":is_in_contact",false,"is_in_contact"),
  lumped_boundary(0,2,0,id+":lumped_boundary",
  std::numeric_limits<Real>::quiet_NaN(),"lumped_boundary"),
  impedance(0,1,0,id+":impedance",std::numeric_limits<Real>::quiet_NaN(),"impedance"),
  contact_surfaces()
{
  AKANTU_DEBUG_IN();

  FEM & boundary_fem = this->model.getFEMBoundary();
  boundary_fem.initShapeFunctions();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNContact::addSurfacePair(const Surface & slave, const Surface & master, UInt surface_normal_dir) {
  AKANTU_DEBUG_IN();

  UInt dim = this->model.getSpatialDimension();
  AKANTU_DEBUG_ASSERT(surface_normal_dir < dim, "Model is of " << dim << " dimensions"
		      << " and cannot have direction " << surface_normal_dir
		      << " for surface normal");

  const Mesh & mesh_ref = this->model.getFEM().getMesh();
  
  const SubBoundary & slave_boundary = mesh_ref.getSubBoundary(slave);
  const SubBoundary & master_boundary = mesh_ref.getSubBoundary(master);
  
  this->contact_surfaces.insert(&slave_boundary);
  this->contact_surfaces.insert(&master_boundary);

  // offset for projection computation
  UInt offset[dim-1];
  for (UInt i=0, j=0; i<dim; ++i) {
    if (surface_normal_dir != i) {
      offset[j] = i;
      ++j;
    }
  }

  // find projected node coordinates
  const Array<Real> & coordinates = mesh_ref.getNodes();
  
  // find slave nodes
  Array<Real> proj_slave_coord(slave_boundary.getNbNodes(),dim-1,0.);
  Array<UInt> slave_nodes(slave_boundary.getNbNodes());
  UInt n(0);
  for(SubBoundary::nodes_const_iterator nodes_it(slave_boundary.nodes_begin()); nodes_it!= slave_boundary.nodes_end(); ++nodes_it) {
    for (UInt d=0; d<dim-1; ++d) {
      proj_slave_coord(n,d) = coordinates(*nodes_it,offset[d]);
      slave_nodes(n)=*nodes_it;
    }
    ++n;
  }
  
  // find master nodes
  Array<Real> proj_master_coord(master_boundary.getNbNodes(),dim-1,0.);
  Array<UInt> master_nodes(master_boundary.getNbNodes());
  n=0;
  for(SubBoundary::nodes_const_iterator nodes_it(master_boundary.nodes_begin()); nodes_it!= master_boundary.nodes_end(); ++nodes_it) {
    for (UInt d=0; d<dim-1; ++d) {
      proj_master_coord(n,d) = coordinates(*nodes_it,offset[d]);
      master_nodes(n)=*nodes_it;
    }
    ++n;
  }

  // find minimum distance between slave nodes to define tolerance
  Real min_dist = std::numeric_limits<Real>::max();
  for (UInt i=0; i<proj_slave_coord.getSize(); ++i) {
    for (UInt j=i+1; j<proj_slave_coord.getSize(); ++j) {
      Real dist = 0.;
      for (UInt d=0; d<dim-1; ++d) {
	      dist += (proj_slave_coord(i,d) - proj_slave_coord(j,d)) 
	      * (proj_slave_coord(i,d) - proj_slave_coord(j,d));
      }
      if (dist < min_dist) {
	      min_dist = dist;
	    }
    }
  }
  min_dist = std::sqrt(min_dist);
  Real local_tol = 0.1*min_dist;

  // find master slave node pairs
  for (UInt i=0; i<proj_slave_coord.getSize(); ++i) {
    for (UInt j=0; j<proj_master_coord.getSize(); ++j) {
      Real dist = 0.;
      for (UInt d=0; d<dim-1; ++d) {
	      dist += (proj_slave_coord(i,d) - proj_master_coord(j,d)) 
	      * (proj_slave_coord(i,d) - proj_master_coord(j,d));
      }
      dist = std::sqrt(dist);
      if (dist < local_tol) { // it is a pair
	      this->addNodePair(slave_nodes(i), master_nodes(j));
	      continue; // found master do not need to search further for this slave
      }
    }
  }

  // synchronize with depending nodes
  updateInternalData();
  syncArrays(_added);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNContact::addNodePair(UInt slave, UInt master) {
  AKANTU_DEBUG_IN();

  UInt dim = this->model.getSpatialDimension();
  
  // add nodes to the node arrays
  this->slaves.push_back(slave);
  this->masters.push_back(master);

  // set contact as false
  this->is_in_contact.push_back(false);

  // before initializing 
  // set contact pressure, impedance, normal, lumped_boundary to Nan
  this->contact_pressure.push_back(std::numeric_limits<Real>::quiet_NaN());
  this->impedance.push_back(std::numeric_limits<Real>::quiet_NaN());
  
  Real nan_two[2];
  nan_two[0] = std::numeric_limits<Real>::quiet_NaN();
  nan_two[1] = std::numeric_limits<Real>::quiet_NaN();
  this->lumped_boundary.push_back(nan_two);
  
  Real nan_normal[dim];
  for (UInt d=0; d<dim; ++d)
    nan_normal[d] = std::numeric_limits<Real>::quiet_NaN();
  this->normals.push_back(nan_normal);

  AKANTU_DEBUG_OUT();  
}

/* -------------------------------------------------------------------------- */
void NTNContact::addNodePairs(Array<UInt> & pairs) {
  AKANTU_DEBUG_IN();
  
  AKANTU_DEBUG_ASSERT(pairs.getNbComponent() == 2, 
		      "Array of node pairs should have nb_component = 2," << 
		      " but has nb_component = " << pairs.getNbComponent());
  UInt nb_pairs = pairs.getSize();
  for (UInt n=0; n<nb_pairs; ++n) {
    this->addNodePair(pairs(n,0), pairs(n,1));
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/*
  This function only works for surface elements with one quad point. For 
  surface elements with more quad points, it computes still, but the result 
  might not be what you are looking for.
 */
void NTNContact::updateNormals() {
  AKANTU_DEBUG_IN();

  UInt dim = this->model.getSpatialDimension();
  const Array<Real> & cur_pos = this->model.getCurrentPosition();
  FEM & boundary_fem = this->model.getFEMBoundary();

  // contact information
  UInt nb_contact_nodes = this->masters.getSize();

  // set normals to zero
  this->normals.clear();

  const Mesh & mesh = this->model.getFEM().getMesh();
  Mesh::type_iterator it = mesh.firstType(dim-1);
  Mesh::type_iterator last = mesh.lastType(dim-1);
  for (; it != last; ++it) {
    // get elements connected to each node
    CSR<UInt> node_to_element;
    MeshUtils::buildNode2ElementsByElementType(mesh, *it, node_to_element);
    // compute the normals
    Array<Real> quad_normals(0,dim);
    boundary_fem.computeNormalsOnControlPoints(cur_pos, quad_normals, *it);

    UInt nb_quad_points = boundary_fem.getNbQuadraturePoints(*it);

    // add up normals to all master nodes
    for (UInt n=0; n<nb_contact_nodes; ++n) {
      UInt master = this->masters(n);
      CSR<UInt>::iterator elem = node_to_element.begin(master);
      // loop over all elements connected to this master node
      for (; elem != node_to_element.end(master); ++elem) {
	UInt e = *elem;
	// loop over all quad points of this element
	for (UInt q=0; q<nb_quad_points; ++q) {
	  // add quad normal to master normal
	  for (UInt d=0; d<dim; ++d) {
	    this->normals(n,d) += quad_normals(e*nb_quad_points + q, d);
	  }
	}
      }
    }
  }

  // normalize normals
  Real * master_normals = this->normals.storage();
  for (UInt n=0; n<nb_contact_nodes; ++n) {
    if (dim==2) 
      Math::normalize2(&(master_normals[n*dim]));
    else if (dim==3)
      Math::normalize3(&(master_normals[n*dim]));
  }
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNContact::registerSyncronizedArray(SyncronizedArrayBase & array) {
  AKANTU_DEBUG_IN();
  
  this->slaves.registerDependingArray(array);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNContact::dumpRestart(const std::string & file_name) const {
  AKANTU_DEBUG_IN();
  
  this->slaves.dumpRestartFile(file_name);
  this->masters.dumpRestartFile(file_name);
  this->normals.dumpRestartFile(file_name);
  this->is_in_contact.dumpRestartFile(file_name);
  this->contact_pressure.dumpRestartFile(file_name);
  this->lumped_boundary.dumpRestartFile(file_name);
  this->impedance.dumpRestartFile(file_name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNContact::readRestart(const std::string & file_name) {
  AKANTU_DEBUG_IN();
  
  this->slaves.readRestartFile(file_name);
  this->masters.readRestartFile(file_name);
  this->normals.readRestartFile(file_name);
  this->is_in_contact.readRestartFile(file_name);
  this->contact_pressure.readRestartFile(file_name);
  this->lumped_boundary.readRestartFile(file_name);
  this->impedance.readRestartFile(file_name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNContact::updateInternalData() {
  AKANTU_DEBUG_IN();
  
  updateNormals();
  updateLumpedBoundary();
  updateImpedance();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNContact::updateImpedance() {
  AKANTU_DEBUG_IN();

  UInt nb_ntn_pairs = getNbContactNodes();
  Real delta_t = this->model.getTimeStep();
  AKANTU_DEBUG_ASSERT(delta_t != NAN, "Time step is NAN. Have you set it already?");

  const Array<Real> & mass = this->model.getMass();

  for (UInt n=0; n<nb_ntn_pairs; ++n) {  
    UInt master = this->masters(n);
    UInt slave = this->slaves(n);

    Real imp = (this->lumped_boundary(n,0) / mass(master)) + (this->lumped_boundary(n,1) / mass(slave));
    imp = 2 / delta_t / imp;
    this->impedance(n) = imp;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNContact::updateLumpedBoundary() {
  AKANTU_DEBUG_IN();

  // set all values in lumped_boundary to zero
  this->lumped_boundary.clear();

  UInt dim = this->model.getSpatialDimension(); 
  UInt nb_contact_nodes = getNbContactNodes();

  const FEM & boundary_fem = this->model.getFEMBoundary();

  const Mesh & mesh = this->model.getFEM().getMesh();
  Mesh::type_iterator it = mesh.firstType(dim-1);
  Mesh::type_iterator last = mesh.lastType(dim-1);
  for (; it != last; ++it) {
    // get elements connected to each node
    CSR<UInt> node_to_element;
    MeshUtils::buildNode2ElementsByElementType(mesh, *it, node_to_element);

    UInt nb_nodes_per_element = mesh.getNbNodesPerElement(*it);
    const Array<UInt> & connectivity = mesh.getConnectivity(*it);

    // get shapes and compute integral 
    const Array<Real> & shapes = boundary_fem.getShapes(*it);
    Array<Real> area(mesh.getNbElement(*it),nb_nodes_per_element);
    boundary_fem.integrate(shapes,area,nb_nodes_per_element,*it);

    // get surface id information
//    const Array<UInt> & surface_id = mesh.getSurfaceID(*it);
//    std::set<UInt>::iterator pos;
//    std::set<UInt>::iterator end = this->contact_surfaces.end();

    // loop over contact nodes
    for (UInt i=0; i<2*nb_contact_nodes; ++i) {
      // first round on masters, second round on slaves
      UInt node = std::numeric_limits<UInt>::quiet_NaN();
      UInt n = 0;
      UInt index = 0;
      if (i<nb_contact_nodes) {
	      n = i;
	      node = this->masters(n);
	      index = 0;
      }
      else {
	      n = i-nb_contact_nodes;
	      node = this->slaves(n);
	      index = 1;
      }

      CSR<UInt>::iterator elem = node_to_element.begin(node);
      // loop over all elements connected to this node
      for (; elem != node_to_element.end(node); ++elem) {
	      UInt e = *elem;

	      // if element is not at interface continue
//	      pos = this->contact_surfaces.find(surface_id(e));
//	      if (pos == end)
//	        continue;

	      // loop over all points of this element
	      for (UInt q=0; q<nb_nodes_per_element; ++q) {
	        if (connectivity(e,q) == node) {
	          this->lumped_boundary(n,index) += area(e,q);
	        }
	      }
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNContact::computeContactPressure() {
  AKANTU_DEBUG_IN();

  UInt dim = this->model.getSpatialDimension();
  Real delta_t = this->model.getTimeStep();
  UInt nb_ntn_pairs = getNbContactNodes();

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
  /*
  for (UInt i=0; i<nb_ntn_pairs; ++i) {
    std::cout << "Master residual = " << residual(masters(i),0) << ", " << residual(masters(i),1) << std::endl;
    std::cout << "Master acceleration = " << acceleration(masters(i),0) << ", " << acceleration(masters(i),1) << std::endl;
    std::cout << "Slave residual = " << residual(slaves(i),0) << ", " << residual(slaves(i),1) << std::endl;
    std::cout << "Slave acceleration = " << acceleration(slaves(i),0) << ", " << acceleration(slaves(i),1) << std::endl;
  }
  */

  // compute relative normal fields of displacement, velocity and acceleration 
  Array<Real> r_disp(0,1);
  Array<Real> r_velo(0,1);
  Array<Real> r_acce(0,1);
  Array<Real> r_old_acce(0,1);
  computeRelativeNormalField(this->model.getCurrentPosition(), r_disp);
  computeRelativeNormalField(this->model.getVelocity(),        r_velo);
  computeRelativeNormalField(acceleration,                     r_acce);
  computeRelativeNormalField(this->model.getAcceleration(),    r_old_acce);

  AKANTU_DEBUG_ASSERT(r_disp.getSize() == nb_ntn_pairs, 
		      "computeRelativeNormalField does not give back arrays " 
		      << "size == nb_ntn_pairs. nb_ntn_pairs = " << nb_ntn_pairs
		      << " | array size = " << r_disp.getSize());

  // compute gap array for all nodes
  Array<Real> gap(nb_ntn_pairs, 1);
  Real * gap_p = gap.storage();
  Real * r_disp_p = r_disp.storage();
  Real * r_velo_p = r_velo.storage();
  Real * r_acce_p = r_acce.storage();
  Real * r_old_acce_p = r_old_acce.storage();
  for (UInt i=0; i<nb_ntn_pairs; ++i) {
    //std::cout << "gap elements: " << *r_disp_p << " " << *r_velo_p << " " << *r_acce_p << std::endl;
    *gap_p = *r_disp_p + delta_t * *r_velo_p + delta_t * delta_t * *r_acce_p - 0.5 * delta_t * delta_t * *r_old_acce_p;
    // increment pointers
    gap_p++;
    r_disp_p++;
    r_velo_p++;
    r_acce_p++;
    r_old_acce_p++;
  }

  // check if gap is negative -> is in contact
  for (UInt n=0; n<nb_ntn_pairs; ++n) {
    if (gap(n) < 0.) {
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
void NTNContact::applyContactPressure() {
  AKANTU_DEBUG_IN();

  UInt nb_ntn_pairs = getNbContactNodes();
  UInt dim = this->model.getSpatialDimension();

  Array<Real> & residual = this->model.getResidual();

  for (UInt n=0; n<nb_ntn_pairs; ++n) {
    UInt master = this->masters(n);
    UInt slave = this->slaves(n);
    
    for (UInt d=0; d<dim; ++d) {
      //if (n==0 && d==1) std::cout << std::setprecision(14) << "mas_res=" << residual(master,d) << " slv_res=" << residual(slave,d) << " " << std::abs(residual(slave,d))-std::abs(residual(master,d)) << std::endl;
      residual(master,d) += this->lumped_boundary(n,0) * this->contact_pressure(n,d);
      residual(slave, d) -= this->lumped_boundary(n,1) * this->contact_pressure(n,d);
      
      /*
      if (n==0 && d==1) {
	std::cout << std::setprecision(14) << "mas_b=" << lumped_boundary(n,0) << " skv_b=" << lumped_boundary(n,1) << std::endl;
	std::cout << std::setprecision(14) << "mas+=" << lumped_boundary(n,0)*contact_pressure(n,d) << " slv+=" << lumped_boundary(n,1)*contact_pressure(n,d) << std::endl;
	std::cout << std::setprecision(14) << "mas_res=" << residual(master,d) << " slv_res=" << residual(slave,d) << std::endl;
	}*/
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNContact::computeRelativeTangentialField(const Array<Real> & field,
						Array<Real> & rel_tang_field) const {
  AKANTU_DEBUG_IN();
  
  // resize arrays to zero
  rel_tang_field.resize(0);
  
  UInt dim = this->model.getSpatialDimension();
  Real * field_p = field.storage();
  Real * normals_p = this->normals.storage();

  UInt nb_contact_nodes = this->slaves.getSize();
  for (UInt n=0; n<nb_contact_nodes; ++n) {
    // nodes
    UInt slave = this->slaves(n);
    UInt master = this->masters(n);

    // compute relative field to master
    Real rel_array[dim];
    for (UInt d=0; d<dim; ++d) {
      rel_array[d] = field_p[slave*dim + d] - field_p[master*dim + d];
    }

    // compute dot product with normal of master
    Real dot_prod = Math::vectorDot(&(normals_p[n*dim]), rel_array, dim);
    
    // compute the tangential projection of the relative field to the master
    for (UInt d=0; d<dim; ++d)
      rel_array[d] -= dot_prod * normals_p[n*dim +d];
    
    rel_tang_field.push_back(rel_array);
  }
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNContact::computeRelativeNormalField(const Array<Real> & field,
					    Array<Real> & rel_normal_field) const {
  AKANTU_DEBUG_IN();
  
  // resize arrays to zero
  rel_normal_field.resize(0);
  
  UInt dim = this->model.getSpatialDimension();
  Real * field_p = field.storage();
  Real * normals_p = this->normals.storage();

  UInt nb_contact_nodes = this->getNbContactNodes();
  for (UInt n=0; n<nb_contact_nodes; ++n) {
    // nodes
    UInt slave = this->slaves(n);
    UInt master = this->masters(n);

    // compute relative field to master
    Real rel_array[dim];
    for (UInt d=0; d<dim; ++d) {
      rel_array[d] = field_p[slave*dim + d] - field_p[master*dim + d];
    }

    // compute dot product with normal of master
    Real dot_prod = Math::vectorDot(&(normals_p[n*dim]), rel_array, dim);
        
    rel_normal_field.push_back(dot_prod);
  }
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Int NTNContact::getNodeIndex(UInt node) const {
  AKANTU_DEBUG_IN();

  Int slave_i  = this->slaves.find(node);
  Int master_i = this->masters.find(node);

  AKANTU_DEBUG_OUT();
  return std::max(slave_i,master_i);
}

/* -------------------------------------------------------------------------- */
void NTNContact::printself(std::ostream & stream, int indent) const {
  AKANTU_DEBUG_IN();
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);
  
  stream << space << "NTNContact [" << std::endl;
  stream << space << " + id            : " << id << std::endl;
  stream << space << " + slaves        : " << std::endl;
  this->slaves.printself(stream, indent + 2);
  stream << space << " + masters       : " << std::endl;
  this->masters.printself(stream, indent + 2);
  stream << space << " + normals       : " << std::endl;
  this->normals.printself(stream, indent + 2);
  stream << space << " + contact_pressure : " << std::endl;
  this->contact_pressure.printself(stream, indent + 2);

  stream << space << "]" << std::endl;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNContact::syncArrays(SyncChoice sync_choice) {
  AKANTU_DEBUG_IN();
  
  this->slaves.syncElements(sync_choice);
  this->masters.syncElements(sync_choice);
  this->normals.syncElements(sync_choice);
  this->is_in_contact.syncElements(sync_choice);
  this->contact_pressure.syncElements(sync_choice);
  this->lumped_boundary.syncElements(sync_choice);
  this->impedance.syncElements(sync_choice);
  
  AKANTU_DEBUG_OUT();
}

__END_SIMTOOLS__
