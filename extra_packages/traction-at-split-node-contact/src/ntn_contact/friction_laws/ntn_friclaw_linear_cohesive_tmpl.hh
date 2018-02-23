/**
 * @file   ntn_friclaw_linear_cohesive_tmpl.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 *
 * @brief  implementation of linear cohesive law
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
//#include "dumper_text.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <class Regularisation>
NTNFricLawLinearCohesive<Regularisation>::NTNFricLawLinearCohesive(
    NTNBaseContact * contact, const FrictionID & id, const MemoryID & memory_id)
    : Regularisation(contact, id, memory_id),
      G_c(0, 1, 0., id + ":G_c", 0., "G_c"),
      tau_c(0, 1, 0., id + ":tau_c", 0., "tau_c"),
      tau_r(0, 1, 0., id + ":tau_r", 0., "tau_r") {
  AKANTU_DEBUG_IN();

  Regularisation::registerSynchronizedArray(this->G_c);
  Regularisation::registerSynchronizedArray(this->tau_c);
  Regularisation::registerSynchronizedArray(this->tau_r);

  this->registerParam("G_c", this->G_c, _pat_parsmod, "fracture energy");
  this->registerParam("tau_c", this->tau_c, _pat_parsmod,
                      "peak shear strength");
  this->registerParam("tau_r", this->tau_r, _pat_parsmod,
                      "residual shear strength");

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawLinearCohesive<Regularisation>::computeFrictionalStrength() {
  AKANTU_DEBUG_IN();

  // get arrays
  const SynchronizedArray<bool> & is_in_contact =
      this->internalGetIsInContact();
  // const SynchronizedArray<Real> & slip  = this->internalGetSlip();
  const SynchronizedArray<Real> & slip = this->internalGetCumulativeSlip();

  // array to fill
  SynchronizedArray<Real> & strength = this->internalGetFrictionalStrength();

  UInt nb_contact_nodes = this->contact->getNbContactNodes();
  for (UInt n = 0; n < nb_contact_nodes; ++n) {
    // node pair is NOT in contact
    if (!is_in_contact(n))
      strength(n) = 0.;

    // node pair is in contact
    else {
      if (this->G_c(n) == 0.) {
        // strength(n) = 0.;
        strength(n) = this->tau_r(n);
      } else {
        Real slope = (this->tau_c(n) - this->tau_r(n)) *
                     (this->tau_c(n) - this->tau_r(n)) / (2 * this->G_c(n));
        // no strength < tau_r
        strength(n) =
            std::max(this->tau_c(n) - slope * slip(n), this->tau_r(n));
        // strength(n) = std::max(this->tau_c(n) - slope * slip(n), 0.); // no
        // negative strength
      }
    }
  }

  Regularisation::computeFrictionalStrength();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawLinearCohesive<Regularisation>::registerSynchronizedArray(
    SynchronizedArrayBase & array) {
  AKANTU_DEBUG_IN();

  this->G_c.registerDependingArray(array);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawLinearCohesive<Regularisation>::dumpRestart(
    const std::string & file_name) const {
  AKANTU_DEBUG_IN();

  this->G_c.dumpRestartFile(file_name);
  this->tau_c.dumpRestartFile(file_name);
  this->tau_r.dumpRestartFile(file_name);

  Regularisation::dumpRestart(file_name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawLinearCohesive<Regularisation>::readRestart(
    const std::string & file_name) {
  AKANTU_DEBUG_IN();

  this->G_c.readRestartFile(file_name);
  this->tau_c.readRestartFile(file_name);
  this->tau_r.readRestartFile(file_name);

  Regularisation::readRestart(file_name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawLinearCohesive<Regularisation>::printself(std::ostream & stream,
                                                         int indent) const {
  AKANTU_DEBUG_IN();
  std::string space;
  for (Int i = 0; i < indent; i++, space += AKANTU_INDENT)
    ;

  stream << space << "NTNFricLawLinearCohesive [" << std::endl;
  Regularisation::printself(stream, ++indent);
  stream << space << "]" << std::endl;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawLinearCohesive<Regularisation>::addDumpFieldToDumper(
    const std::string & dumper_name, const std::string & field_id) {
  AKANTU_DEBUG_IN();

#ifdef AKANTU_USE_IOHELPER
  //  const SynchronizedArray<UInt> * nodal_filter =
  //  &(this->contact->getSlaves());

  if (field_id == "G_c") {
    this->internalAddDumpFieldToDumper(
        dumper_name, field_id,
        new dumper::NodalField<Real>(this->G_c.getArray()));
  } else if (field_id == "tau_c") {
    this->internalAddDumpFieldToDumper(
        dumper_name, field_id,
        new dumper::NodalField<Real>(this->tau_c.getArray()));
  } else if (field_id == "tau_r") {
    this->internalAddDumpFieldToDumper(
        dumper_name, field_id,
        new dumper::NodalField<Real>(this->tau_r.getArray()));
  } else {
    Regularisation::addDumpFieldToDumper(dumper_name, field_id);
  }

#endif

  AKANTU_DEBUG_OUT();
}

} // namespace akantu
