/**
 * @file   ntn_base_contact.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Tue Dec 02 2014
 * @date last modification: Fri Feb 23 2018
 *
 * @brief  implementation of ntn base contact
 *
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "ntn_base_contact.hh"
#include "dof_manager_default.hh"
#include "dumpable_inline_impl.hh"
#include "dumper_nodal_field.hh"
#include "dumper_text.hh"
#include "element_synchronizer.hh"
#include "mesh_utils.hh"
#include "non_linear_solver_lumped.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
// NTNContactSynchElementFilter::NTNContactSynchElementFilter(
//     NTNBaseContact & contact)
//     : contact(contact),
//       connectivity(contact.getModel().getMesh().getConnectivities()) {
//   AKANTU_DEBUG_IN();

//   AKANTU_DEBUG_OUT();
// }

/* -------------------------------------------------------------------------- */
// bool NTNContactSynchElementFilter::operator()(const Element & e) {
//   AKANTU_DEBUG_IN();

//   ElementType type = e.type;
//   UInt element = e.element;
//   GhostType ghost_type = e.ghost_type;

//   // loop over all nodes of this element
//   bool need_element = false;
//   UInt nb_nodes = Mesh::getNbNodesPerElement(type);
//   for (UInt n = 0; n < nb_nodes; ++n) {
//     UInt nn = this->connectivity(type, ghost_type)(element, n);

//     // if one nodes is in this contact, we need this element
//     if (this->contact.getNodeIndex(nn) >= 0) {
//       need_element = true;
//       break;
//     }
//   }

//   AKANTU_DEBUG_OUT();
//   return need_element;
// }

/* -------------------------------------------------------------------------- */
NTNBaseContact::NTNBaseContact(SolidMechanicsModel & model, const ID & id,
                               const MemoryID & memory_id)
    : Memory(id, memory_id), Dumpable(), model(model),
      slaves(0, 1, 0, id + ":slaves", std::numeric_limits<UInt>::quiet_NaN(),
             "slaves"),
      normals(0, model.getSpatialDimension(), 0, id + ":normals",
              std::numeric_limits<Real>::quiet_NaN(), "normals"),
      contact_pressure(
          0, model.getSpatialDimension(), 0, id + ":contact_pressure",
          std::numeric_limits<Real>::quiet_NaN(), "contact_pressure"),
      is_in_contact(0, 1, false, id + ":is_in_contact", false, "is_in_contact"),
      lumped_boundary_slaves(0, 1, 0, id + ":lumped_boundary_slaves",
                             std::numeric_limits<Real>::quiet_NaN(),
                             "lumped_boundary_slaves"),
      impedance(0, 1, 0, id + ":impedance",
                std::numeric_limits<Real>::quiet_NaN(), "impedance"),
      contact_surfaces(), slave_elements("slave_elements", id, memory_id),
      node_to_elements() {
  AKANTU_DEBUG_IN();

  auto & boundary_fem = this->model.getFEEngineBoundary();
  for (auto && ghost_type : ghost_types) {
    boundary_fem.initShapeFunctions(ghost_type);
  }

  auto & mesh = this->model.getMesh();
  auto spatial_dimension = this->model.getSpatialDimension();

  this->slave_elements.initialize(mesh,
                                  _spatial_dimension = spatial_dimension - 1);

  MeshUtils::buildNode2Elements(mesh, this->node_to_elements,
                                spatial_dimension - 1);

  this->registerDumper<DumperText>("text_all", id, true);
  this->addDumpFilteredMesh(mesh, slave_elements, slaves.getArray(),
                            spatial_dimension - 1, _not_ghost, _ek_regular);

  // parallelisation
  this->synch_registry = std::make_unique<SynchronizerRegistry>();
  this->synch_registry->registerDataAccessor(*this);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
NTNBaseContact::~NTNBaseContact() = default;

/* -------------------------------------------------------------------------- */
void NTNBaseContact::initParallel() {
  AKANTU_DEBUG_IN();

  this->synchronizer = std::make_unique<ElementSynchronizer>(
      this->model.getMesh().getElementSynchronizer());

  this->synchronizer->filterScheme([&](auto && element) {
    // loop over all nodes of this element
    Vector<UInt> conn = const_cast<const Mesh &>(this->model.getMesh())
                            .getConnectivity(element);
    for (auto & node : conn) {
      // if one nodes is in this contact, we need this element
      if (this->getNodeIndex(node) >= 0) {
        return true;
      }
    }
    return false;
  });

  this->synch_registry->registerSynchronizer(*(this->synchronizer),
                                             SynchronizationTag::_cf_nodal);
  this->synch_registry->registerSynchronizer(*(this->synchronizer),
                                             SynchronizationTag::_cf_incr);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNBaseContact::findBoundaryElements(
    const Array<UInt> & interface_nodes, ElementTypeMapArray<UInt> & elements) {
  AKANTU_DEBUG_IN();

  // add connected boundary elements that have all nodes on this contact
  for (const auto & node : interface_nodes) {
    for (const auto & element : this->node_to_elements.getRow(node)) {
      Vector<UInt> conn = const_cast<const Mesh &>(this->model.getMesh())
                              .getConnectivity(element);
      auto nb_nodes = conn.size();
      decltype(nb_nodes) nb_found_nodes = 0;

      for (auto & nn : conn) {
        if (interface_nodes.find(nn) != UInt(-1)) {
          nb_found_nodes++;
        } else {
          break;
        }
      }

      // this is an element between all contact nodes
      // and is not already in the elements
      if ((nb_found_nodes == nb_nodes) &&
          (elements(element.type, element.ghost_type).find(element.element) ==
           UInt(-1))) {
        elements(element.type, element.ghost_type).push_back(element.element);
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNBaseContact::addSplitNode(UInt node, UInt) {
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
  this->lumped_boundary_slaves.push_back(
      std::numeric_limits<Real>::quiet_NaN());

  Vector<Real> nan_normal(dim, std::numeric_limits<Real>::quiet_NaN());
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
    bool is_pbc_slave_node = mesh.isPeriodicSlave(this->slaves(n));
    if (is_local_node && !is_pbc_slave_node && this->is_in_contact(n)) {
      nb_contact++;
    }
  }

  mesh.getCommunicator().allReduce(nb_contact, SynchronizerOperation::_sum);

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
void NTNBaseContact::internalUpdateLumpedBoundary(
    const Array<UInt> & nodes, const ElementTypeMapArray<UInt> & elements,
    SynchronizedArray<Real> & boundary) {
  AKANTU_DEBUG_IN();

  // set all values in lumped_boundary to zero
  boundary.zero();

  UInt dim = this->model.getSpatialDimension();
  //  UInt nb_contact_nodes = getNbContactNodes();

  const FEEngine & boundary_fem = this->model.getFEEngineBoundary();

  const Mesh & mesh = this->model.getMesh();

  for (auto ghost_type : ghost_types) {
    for (auto & type : mesh.elementTypes(dim - 1, ghost_type)) {
      UInt nb_elements = mesh.getNbElement(type, ghost_type);
      UInt nb_nodes_per_element = mesh.getNbNodesPerElement(type);
      const Array<UInt> & connectivity = mesh.getConnectivity(type, ghost_type);

      // get shapes and compute integral
      const Array<Real> & shapes = boundary_fem.getShapes(type, ghost_type);
      Array<Real> area(nb_elements, nb_nodes_per_element);
      boundary_fem.integrate(shapes, area, nb_nodes_per_element, type,
                             ghost_type);

      if (this->contact_surfaces.size() == 0) {
        AKANTU_DEBUG_WARNING(
            "No surfaces in ntn base contact."
            << " You have to define the lumped boundary by yourself.");
      }

      Array<UInt>::const_iterator<UInt> elem_it =
          (elements)(type, ghost_type).begin();
      Array<UInt>::const_iterator<UInt> elem_it_end =
          (elements)(type, ghost_type).end();
      // loop over contact nodes
      for (; elem_it != elem_it_end; ++elem_it) {
        for (UInt q = 0; q < nb_nodes_per_element; ++q) {
          UInt node = connectivity(*elem_it, q);
          UInt node_index = nodes.find(node);
          AKANTU_DEBUG_ASSERT(node_index != UInt(-1), "Could not find node "
                                                          << node
                                                          << " in the array!");
          Real area_to_add = area(*elem_it, q);
          boundary(node_index) += area_to_add;
        }
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNBaseContact::computeAcceleration(Array<Real> & acceleration) const {
  auto && dof_manager =
      dynamic_cast<DOFManagerDefault &>(model.getDOFManager());
  const auto & b = dof_manager.getResidual();
  acceleration.resize(b.size());
  const auto & blocked_dofs = dof_manager.getGlobalBlockedDOFs();
  const auto & A = dof_manager.getLumpedMatrix("M");

  Array<bool> blocked_dofs_bool(blocked_dofs.size());
  for (auto && data : zip(blocked_dofs, blocked_dofs_bool)) {
    std::get<1>(data) = std::get<0>(data);
  }

  // pre-compute the acceleration
  // (not increment acceleration, because residual is still Kf)
  NonLinearSolverLumped::solveLumped(A, acceleration, b, this->model.getF_M2A(),
                                     blocked_dofs_bool);
}

/* -------------------------------------------------------------------------- */
void NTNBaseContact::computeContactPressure() {
  AKANTU_DEBUG_IN();

  UInt dim = this->model.getSpatialDimension();
  Real delta_t = this->model.getTimeStep();
  UInt nb_contact_nodes = getNbContactNodes();

  AKANTU_DEBUG_ASSERT(delta_t > 0.,
                      "Cannot compute contact pressure if no time step is set");

  // synchronize data
  this->synch_registry->synchronize(SynchronizationTag::_cf_nodal);

  Array<Real> acceleration(0, dim);
  this->computeAcceleration(acceleration);

  // compute relative normal fields of displacement, velocity and acceleration
  Array<Real> r_disp(0, 1);
  Array<Real> r_velo(0, 1);
  Array<Real> r_acce(0, 1);
  Array<Real> r_old_acce(0, 1);
  computeNormalGap(r_disp);
  //  computeRelativeNormalField(this->model.getCurrentPosition(), r_disp);
  computeRelativeNormalField(this->model.getVelocity(), r_velo);
  computeRelativeNormalField(acceleration, r_acce);
  computeRelativeNormalField(this->model.getAcceleration(), r_old_acce);

  AKANTU_DEBUG_ASSERT(r_disp.size() == nb_contact_nodes,
                      "computeRelativeNormalField does not give back arrays "
                          << "size == nb_contact_nodes. nb_contact_nodes = "
                          << nb_contact_nodes
                          << " | array size = " << r_disp.size());

  // compute gap array for all nodes
  Array<Real> gap(nb_contact_nodes, 1);
  Real * gap_p = gap.storage();
  Real * r_disp_p = r_disp.storage();
  Real * r_velo_p = r_velo.storage();
  Real * r_acce_p = r_acce.storage();
  Real * r_old_acce_p = r_old_acce.storage();
  for (UInt i = 0; i < nb_contact_nodes; ++i) {
    *gap_p = *r_disp_p + delta_t * *r_velo_p + delta_t * delta_t * *r_acce_p -
             0.5 * delta_t * delta_t * *r_old_acce_p;
    // increment pointers
    gap_p++;
    r_disp_p++;
    r_velo_p++;
    r_acce_p++;
    r_old_acce_p++;
  }

  // check if gap is negative -> is in contact
  for (UInt n = 0; n < nb_contact_nodes; ++n) {
    if (gap(n) <= 0.) {
      for (UInt d = 0; d < dim; ++d) {
        this->contact_pressure(n, d) =
            this->impedance(n) * gap(n) / (2 * delta_t) * this->normals(n, d);
      }
      this->is_in_contact(n) = true;
    } else {
      for (UInt d = 0; d < dim; ++d)
        this->contact_pressure(n, d) = 0.;
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

  Array<Real> & residual = this->model.getInternalForce();

  for (UInt n = 0; n < nb_contact_nodes; ++n) {
    UInt slave = this->slaves(n);

    for (UInt d = 0; d < dim; ++d) {
      // residual(master,d) += this->lumped_boundary(n,0) *
      // this->contact_pressure(n,d);
      residual(slave, d) -=
          this->lumped_boundary_slaves(n) * this->contact_pressure(n, d);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Int NTNBaseContact::getNodeIndex(UInt node) const {
  return this->slaves.find(node);
}

/* -------------------------------------------------------------------------- */
void NTNBaseContact::printself(std::ostream & stream, int indent) const {
  AKANTU_DEBUG_IN();
  std::string space;
  for (Int i = 0; i < indent; i++, space += AKANTU_INDENT)
    ;

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

#define ADD_FIELD(field_id, field, type)                                       \
  internalAddDumpFieldToDumper(                                                \
      dumper_name, field_id,                                                   \
      std::make_unique<                                                        \
          dumpers::NodalField<type, true, Array<type>, Array<UInt>>>(           \
          field, 0, 0, &nodal_filter))

  if (field_id == "displacement") {
    ADD_FIELD(field_id, this->model.getDisplacement(), Real);
  } else if (field_id == "mass") {
    ADD_FIELD(field_id, this->model.getMass(), Real);
  } else if (field_id == "velocity") {
    ADD_FIELD(field_id, this->model.getVelocity(), Real);
  } else if (field_id == "acceleration") {
    ADD_FIELD(field_id, this->model.getAcceleration(), Real);
  } else if (field_id == "external_force") {
    ADD_FIELD(field_id, this->model.getExternalForce(), Real);
  } else if (field_id == "internal_force") {
    ADD_FIELD(field_id, this->model.getInternalForce(), Real);
  } else if (field_id == "blocked_dofs") {
    ADD_FIELD(field_id, this->model.getBlockedDOFs(), bool);
  } else if (field_id == "increment") {
    ADD_FIELD(field_id, this->model.getIncrement(), Real);
  } else if (field_id == "normal") {
    internalAddDumpFieldToDumper(
        dumper_name, field_id,
        std::make_unique<dumpers::NodalField<Real>>(this->normals.getArray()));
  } else if (field_id == "contact_pressure") {
    internalAddDumpFieldToDumper(dumper_name, field_id,
                                 std::make_unique<dumpers::NodalField<Real>>(
                                     this->contact_pressure.getArray()));
  } else if (field_id == "is_in_contact") {
    internalAddDumpFieldToDumper(dumper_name, field_id,
                                 std::make_unique<dumpers::NodalField<bool>>(
                                     this->is_in_contact.getArray()));
  } else if (field_id == "lumped_boundary_slave") {
    internalAddDumpFieldToDumper(dumper_name, field_id,
                                 std::make_unique<dumpers::NodalField<Real>>(
                                     this->lumped_boundary_slaves.getArray()));
  } else if (field_id == "impedance") {
    internalAddDumpFieldToDumper(
        dumper_name, field_id,
        std::make_unique<dumpers::NodalField<Real>>(this->impedance.getArray()));
  } else {
    std::cerr << "Could not add field '" << field_id
              << "' to the dumper. Just ignored it." << std::endl;
  }

#undef ADD_FIELD
#endif

  AKANTU_DEBUG_OUT();
}

} // namespace akantu
