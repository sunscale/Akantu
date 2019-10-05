/**
 * @file   ntn_friclaw_coulomb_tmpl.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Tue Dec 02 2014
 * @date last modification: Fri Feb 23 2018
 *
 * @brief  implementation of coulomb friction
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
#include "dumper_nodal_field.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <class Regularisation>
NTNFricLawCoulomb<Regularisation>::NTNFricLawCoulomb(NTNBaseContact & contact,
                                                     const ID & id,
                                                     const MemoryID & memory_id)
    : Regularisation(contact, id, memory_id),
      mu(0, 1, 0., id + ":mu", 0., "mu") {
  AKANTU_DEBUG_IN();

  Regularisation::registerSynchronizedArray(this->mu);

  this->registerParam("mu", this->mu, _pat_parsmod, "friction coefficient");

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawCoulomb<Regularisation>::computeFrictionalStrength() {
  AKANTU_DEBUG_IN();

  // get contact arrays
  const SynchronizedArray<bool> & is_in_contact =
      this->internalGetIsInContact();
  const SynchronizedArray<Real> & pressure = this->internalGetContactPressure();

  // array to fill
  SynchronizedArray<Real> & strength = this->internalGetFrictionalStrength();

  UInt nb_contact_nodes = this->contact.getNbContactNodes();
  for (UInt n = 0; n < nb_contact_nodes; ++n) {
    // node pair is NOT in contact
    if (!is_in_contact(n))
      strength(n) = 0.;

    // node pair is in contact
    else {
      // compute frictional strength
      strength(n) = this->mu(n) * pressure(n);
    }
  }

  Regularisation::computeFrictionalStrength();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawCoulomb<Regularisation>::registerSynchronizedArray(
    SynchronizedArrayBase & array) {
  AKANTU_DEBUG_IN();

  this->mu.registerDependingArray(array);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawCoulomb<Regularisation>::dumpRestart(
    const std::string & file_name) const {
  AKANTU_DEBUG_IN();

  this->mu.dumpRestartFile(file_name);

  Regularisation::dumpRestart(file_name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawCoulomb<Regularisation>::readRestart(
    const std::string & file_name) {
  AKANTU_DEBUG_IN();

  this->mu.readRestartFile(file_name);

  Regularisation::readRestart(file_name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawCoulomb<Regularisation>::printself(std::ostream & stream,
                                                  int indent) const {
  AKANTU_DEBUG_IN();
  std::string space;
  for (Int i = 0; i < indent; i++, space += AKANTU_INDENT)
    ;

  stream << space << "NTNFricLawCoulomb [" << std::endl;
  Regularisation::printself(stream, ++indent);
  stream << space << "]" << std::endl;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawCoulomb<Regularisation>::addDumpFieldToDumper(
    const std::string & dumper_name, const std::string & field_id) {
  AKANTU_DEBUG_IN();

#ifdef AKANTU_USE_IOHELPER
  //  const SynchronizedArray<UInt> * nodal_filter =
  //  &(this->contact.getSlaves());

  if (field_id == "mu") {
    this->internalAddDumpFieldToDumper(
        dumper_name, field_id,
        std::make_unique<dumper::NodalField<Real>>(this->mu.getArray()));
  }
  /*
  else if (field_id == "frictional_contact_pressure") {
    this->internalAddDumpFieldToDumper(dumper_name,
                       field_id,
                       new
  DumperIOHelper::NodalField<Real>(this->frictional_contact_pressure.getArray()));
  }
  */
  else {
    Regularisation::addDumpFieldToDumper(dumper_name, field_id);
  }

#endif

  AKANTU_DEBUG_OUT();
}

} // namespace akantu
