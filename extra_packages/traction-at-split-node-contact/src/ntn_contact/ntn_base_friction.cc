/**
 * @file   ntn_base_friction.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Tue Dec 02 2014
 * @date last modification: Fri Feb 23 2018
 *
 * @brief  implementation of ntn base friction
 *
 * @section LICENSE
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
// simtools
#include "ntn_base_friction.hh"
#include "dof_manager_default.hh"
#include "dumper_nodal_field.hh"
#include "dumper_text.hh"
#include "non_linear_solver_lumped.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
NTNBaseFriction::NTNBaseFriction(NTNBaseContact & contact, const ID & id,
                                 const MemoryID & memory_id)
    : Memory(id, memory_id), Parsable(ParserType::_friction, id), Dumpable(),
      contact(contact),
      is_sticking(0, 1, true, id + ":is_sticking", true, "is_sticking"),
      frictional_strength(0, 1, 0., id + ":frictional_strength", 0.,
                          "frictional_strength"),
      friction_traction(0, contact.getModel().getSpatialDimension(), 0.,
                        id + ":friction_traction", 0., "friction_traction"),
      slip(0, 1, 0., id + ":slip", 0., "slip"),
      cumulative_slip(0, 1, 0., id + ":cumulative_slip", 0., "cumulative_slip"),
      slip_velocity(0, contact.getModel().getSpatialDimension(), 0.,
                    id + ":slip_velocity", 0., "slip_velocity") {
  AKANTU_DEBUG_IN();

  this->contact.registerSynchronizedArray(this->is_sticking);
  this->contact.registerSynchronizedArray(this->frictional_strength);
  this->contact.registerSynchronizedArray(this->friction_traction);
  this->contact.registerSynchronizedArray(this->slip);
  this->contact.registerSynchronizedArray(this->cumulative_slip);
  this->contact.registerSynchronizedArray(this->slip_velocity);

  this->registerExternalDumper(contact.getDumper(),
                               contact.getDefaultDumperName(), true);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNBaseFriction::updateSlip() {
  AKANTU_DEBUG_IN();

  SolidMechanicsModel & model = this->contact.getModel();
  UInt dim = model.getSpatialDimension();

  // synchronize increment
  this->contact.getSynchronizerRegistry().synchronize(
      SynchronizationTag::_cf_incr);

  Array<Real> rel_tan_incr(0, dim);
  this->contact.computeRelativeTangentialField(model.getIncrement(),
                                               rel_tan_incr);
  Array<Real>::const_iterator<Vector<Real>> it = rel_tan_incr.begin(dim);

  UInt nb_nodes = this->contact.getNbContactNodes();
  for (UInt n = 0; n < nb_nodes; ++n) {
    if (this->is_sticking(n)) {
      this->slip(n) = 0.;
    } else {
      const Vector<Real> & rti = it[n];
      this->slip(n) += rti.norm();
      this->cumulative_slip(n) += rti.norm();
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNBaseFriction::computeFrictionTraction() {
  AKANTU_DEBUG_IN();

  this->computeStickTraction();
  this->computeFrictionalStrength();

  SolidMechanicsModel & model = this->contact.getModel();
  UInt dim = model.getSpatialDimension();

  // get contact arrays
  const SynchronizedArray<bool> & is_in_contact =
      this->contact.getIsInContact();

  Array<Real> & traction =
      const_cast<Array<Real> &>(this->friction_traction.getArray());
  Array<Real>::iterator<Vector<Real>> it_fric_trac = traction.begin(dim);

  this->is_sticking.clear(); // set to not sticking

  UInt nb_contact_nodes = this->contact.getNbContactNodes();
  for (UInt n = 0; n < nb_contact_nodes; ++n) {
    // node pair is in contact
    if (is_in_contact(n)) {
      Vector<Real> fric_trac = it_fric_trac[n];
      // check if it is larger than frictional strength
      Real abs_fric = fric_trac.norm();
      if (abs_fric != 0.) {
        Real alpha = this->frictional_strength(n) / abs_fric;

        // larger -> sliding
        if (alpha < 1.) {
          fric_trac *= alpha;
        } else
          this->is_sticking(n) = true;
      } else {
        // frictional traction is already zero
        this->is_sticking(n) = true;
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNBaseFriction::computeStickTraction() {
  AKANTU_DEBUG_IN();

  SolidMechanicsModel & model = this->contact.getModel();
  UInt dim = model.getSpatialDimension();
  Real delta_t = model.getTimeStep();

  UInt nb_contact_nodes = this->contact.getNbContactNodes();

  // get contact arrays
  const SynchronizedArray<Real> & impedance = this->contact.getImpedance();
  const SynchronizedArray<bool> & is_in_contact =
      this->contact.getIsInContact();

  Array<Real> acceleration(0, dim);
  this->contact.computeAcceleration(acceleration);

  // compute relative normal fields of velocity and acceleration
  Array<Real> r_velo(0, dim);
  Array<Real> r_acce(0, dim);
  Array<Real> r_old_acce(0, dim);
  this->contact.computeRelativeTangentialField(model.getVelocity(), r_velo);
  this->contact.computeRelativeTangentialField(acceleration, r_acce);
  this->contact.computeRelativeTangentialField(model.getAcceleration(),
                                               r_old_acce);

  AKANTU_DEBUG_ASSERT(r_velo.size() == nb_contact_nodes,
                      "computeRelativeNormalField does not give back arrays "
                          << "size == nb_contact_nodes. nb_contact_nodes = "
                          << nb_contact_nodes
                          << " | array size = " << r_velo.size());

  // compute tangential gap_dot array for all nodes
  Array<Real> gap_dot(nb_contact_nodes, dim);
  for (auto && data : zip(make_view(gap_dot), make_view(r_velo),
                          make_view(r_acce), make_view(r_old_acce))) {
    auto & gap_dot = std::get<0>(data);
    auto & r_velo = std::get<1>(data);
    auto & r_acce = std::get<2>(data);
    auto & r_old_acce = std::get<3>(data);

    gap_dot = r_velo + delta_t * r_acce - 1. / 2. * delta_t * r_old_acce;
  }

  // compute friction traction to stop sliding
  Array<Real> & traction =
      const_cast<Array<Real> &>(this->friction_traction.getArray());
  auto it_fric_trac = traction.begin(dim);
  for (UInt n = 0; n < nb_contact_nodes; ++n) {
    Vector<Real> fric_trac = it_fric_trac[n];
    // node pair is NOT in contact
    if (!is_in_contact(n)) {
      fric_trac.clear(); // set to zero
    }

    // node pair is in contact
    else {
      // compute friction traction
      for (UInt d = 0; d < dim; ++d)
        fric_trac(d) = impedance(n) * gap_dot(n, d) / 2.;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNBaseFriction::applyFrictionTraction() {
  AKANTU_DEBUG_IN();

  SolidMechanicsModel & model = this->contact.getModel();
  Array<Real> & residual = model.getInternalForce();
  UInt dim = model.getSpatialDimension();

  const SynchronizedArray<UInt> & slaves = this->contact.getSlaves();
  const SynchronizedArray<Real> & lumped_boundary_slaves =
      this->contact.getLumpedBoundarySlaves();

  UInt nb_contact_nodes = this->contact.getNbContactNodes();
  for (UInt n = 0; n < nb_contact_nodes; ++n) {
    UInt slave = slaves(n);

    for (UInt d = 0; d < dim; ++d) {
      residual(slave, d) -=
          lumped_boundary_slaves(n) * this->friction_traction(n, d);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNBaseFriction::registerSynchronizedArray(SynchronizedArrayBase & array) {
  AKANTU_DEBUG_IN();

  this->frictional_strength.registerDependingArray(array);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNBaseFriction::dumpRestart(const std::string & file_name) const {
  AKANTU_DEBUG_IN();

  this->is_sticking.dumpRestartFile(file_name);
  this->frictional_strength.dumpRestartFile(file_name);
  this->friction_traction.dumpRestartFile(file_name);
  this->slip.dumpRestartFile(file_name);
  this->cumulative_slip.dumpRestartFile(file_name);
  this->slip_velocity.dumpRestartFile(file_name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNBaseFriction::readRestart(const std::string & file_name) {
  AKANTU_DEBUG_IN();

  this->is_sticking.readRestartFile(file_name);
  this->frictional_strength.readRestartFile(file_name);
  this->friction_traction.readRestartFile(file_name);
  this->cumulative_slip.readRestartFile(file_name);
  this->slip_velocity.readRestartFile(file_name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNBaseFriction::setParam(const std::string & name, UInt node,
                               Real value) {
  AKANTU_DEBUG_IN();

  SynchronizedArray<Real> & array =
      this->get(name).get<SynchronizedArray<Real>>();
  Int index = this->contact.getNodeIndex(node);
  if (index < 0) {
    AKANTU_DEBUG_WARNING("Node "
                         << node << " is not a contact node. "
                         << "Therefore, cannot set interface parameter!!");
  } else {
    array(index) = value; // put value
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
UInt NTNBaseFriction::getNbStickingNodes() const {
  AKANTU_DEBUG_IN();

  UInt nb_stick = 0;

  UInt nb_nodes = this->contact.getNbContactNodes();
  const SynchronizedArray<UInt> & nodes = this->contact.getSlaves();
  const SynchronizedArray<bool> & is_in_contact =
      this->contact.getIsInContact();

  const Mesh & mesh = this->contact.getModel().getMesh();

  for (UInt n = 0; n < nb_nodes; ++n) {
    bool is_local_node = mesh.isLocalOrMasterNode(nodes(n));
    bool is_pbc_slave_node = mesh.isPeriodicSlave(nodes(n));
    if (is_local_node && !is_pbc_slave_node && is_in_contact(n) &&
        this->is_sticking(n)) {
      nb_stick++;
    }
  }

  mesh.getCommunicator().allReduce(nb_stick, SynchronizerOperation::_sum);

  AKANTU_DEBUG_OUT();
  return nb_stick;
}

/* -------------------------------------------------------------------------- */
void NTNBaseFriction::printself(std::ostream & stream, int indent) const {
  AKANTU_DEBUG_IN();
  std::string space;
  for (Int i = 0; i < indent; i++, space += AKANTU_INDENT)
    ;

  stream << space << "NTNBaseFriction [" << std::endl;
  Parsable::printself(stream, indent);
  stream << space << "]" << std::endl;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNBaseFriction::addDumpFieldToDumper(const std::string & dumper_name,
                                           const std::string & field_id) {
  AKANTU_DEBUG_IN();

#ifdef AKANTU_USE_IOHELPER
  //  const SynchronizedArray<UInt> * nodal_filter =
  //  &(this->contact.getSlaves());

  if (field_id == "is_sticking") {
    this->internalAddDumpFieldToDumper(
        dumper_name, field_id,
        std::make_unique<dumper::NodalField<bool>>(
            this->is_sticking.getArray()));
  } else if (field_id == "frictional_strength") {
    this->internalAddDumpFieldToDumper(
        dumper_name, field_id,
        std::make_unique<dumper::NodalField<Real>>(
            this->frictional_strength.getArray()));
  } else if (field_id == "friction_traction") {
    this->internalAddDumpFieldToDumper(
        dumper_name, field_id,
        std::make_unique<dumper::NodalField<Real>>(
            this->friction_traction.getArray()));
  } else if (field_id == "slip") {
    this->internalAddDumpFieldToDumper(
        dumper_name, field_id,
        std::make_unique<dumper::NodalField<Real>>(this->slip.getArray()));
  } else if (field_id == "cumulative_slip") {
    this->internalAddDumpFieldToDumper(
        dumper_name, field_id,
        std::make_unique<dumper::NodalField<Real>>(
            this->cumulative_slip.getArray()));
  } else if (field_id == "slip_velocity") {
    this->internalAddDumpFieldToDumper(
        dumper_name, field_id,
        std::make_unique<dumper::NodalField<Real>>(
            this->slip_velocity.getArray()));
  } else {
    this->contact.addDumpFieldToDumper(dumper_name, field_id);
  }

#endif

  AKANTU_DEBUG_OUT();
}

} // namespace akantu
