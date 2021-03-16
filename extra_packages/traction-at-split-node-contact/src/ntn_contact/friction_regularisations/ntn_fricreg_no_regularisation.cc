/**
 * @file   ntn_fricreg_no_regularisation.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Tue Dec 02 2014
 * @date last modification: Fri Feb 23 2018
 *
 * @brief  implementation of no regularisation
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
// simtools
#include "ntn_fricreg_no_regularisation.hh"
#include "dumper_nodal_field.hh"
#include "dumper_text.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
NTNFricRegNoRegularisation::NTNFricRegNoRegularisation(
    NTNBaseContact & contact, const ID & id)
    : NTNBaseFriction(contact, id),
      frictional_contact_pressure(0, 1, 0., id + ":frictional_contact_pressure",
                                  0., "frictional_contact_pressure") {
  AKANTU_DEBUG_IN();

  NTNBaseFriction::registerSynchronizedArray(this->frictional_contact_pressure);

  this->registerParam("frictional_contact_pressure",
                      this->frictional_contact_pressure, _pat_internal,
                      "contact pressure used for friction law");

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
const SynchronizedArray<Real> &
NTNFricRegNoRegularisation::internalGetContactPressure() {
  AKANTU_DEBUG_IN();

  this->computeFrictionalContactPressure();

  AKANTU_DEBUG_OUT();
  return this->frictional_contact_pressure;
}

/* -------------------------------------------------------------------------- */
void NTNFricRegNoRegularisation::computeFrictionalContactPressure() {
  AKANTU_DEBUG_IN();

  SolidMechanicsModel & model = this->contact.getModel();
  UInt dim = model.getSpatialDimension();

  // get contact arrays
  const SynchronizedArray<bool> & is_in_contact =
      this->internalGetIsInContact();
  const Array<Real> & pressure = this->contact.getContactPressure().getArray();
  Array<Real>::const_iterator<Vector<Real>> it = pressure.begin(dim);

  UInt nb_contact_nodes = this->contact.getNbContactNodes();
  for (UInt n = 0; n < nb_contact_nodes; ++n) {
    // node pair is NOT in contact
    if (!is_in_contact(n))
      this->frictional_contact_pressure(n) = 0.;

    // node pair is in contact
    else {
      // compute frictional contact pressure
      const Vector<Real> & pres = it[n];
      this->frictional_contact_pressure(n) = pres.norm();
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFricRegNoRegularisation::registerSynchronizedArray(
    SynchronizedArrayBase & array) {
  AKANTU_DEBUG_IN();

  this->frictional_contact_pressure.registerDependingArray(array);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFricRegNoRegularisation::dumpRestart(
    const std::string & file_name) const {
  AKANTU_DEBUG_IN();

  this->frictional_contact_pressure.dumpRestartFile(file_name);

  NTNBaseFriction::dumpRestart(file_name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFricRegNoRegularisation::readRestart(const std::string & file_name) {
  AKANTU_DEBUG_IN();

  this->frictional_contact_pressure.readRestartFile(file_name);

  NTNBaseFriction::readRestart(file_name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFricRegNoRegularisation::printself(std::ostream & stream,
                                           int indent) const {
  AKANTU_DEBUG_IN();
  std::string space;
  for (Int i = 0; i < indent; i++, space += AKANTU_INDENT)
    ;

  stream << space << "NTNFricRegNoRegularisation [" << std::endl;
  NTNBaseFriction::printself(stream, ++indent);
  stream << space << "]" << std::endl;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFricRegNoRegularisation::addDumpFieldToDumper(
    const std::string & dumper_name, const std::string & field_id) {
  AKANTU_DEBUG_IN();

#ifdef AKANTU_USE_IOHELPER
  //  const SynchronizedArray<UInt> * nodal_filter =
  //  &(this->contact.getSlaves());

  if (field_id == "frictional_contact_pressure") {
    this->internalAddDumpFieldToDumper(
        dumper_name, field_id,
        std::make_unique<dumpers::NodalField<Real>>(
            this->frictional_contact_pressure.getArray()));
  } else {
    NTNBaseFriction::addDumpFieldToDumper(dumper_name, field_id);
  }

#endif

  AKANTU_DEBUG_OUT();
}

} // namespace akantu
