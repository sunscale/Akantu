/**
 * @file   ntn_friclaw_linear_slip_weakening_tmpl.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 *
 * @brief  implementation of linear slip weakening
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#include "dumper_text.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template <class Regularisation>
NTNFricLawLinearSlipWeakening<Regularisation>::NTNFricLawLinearSlipWeakening(NTNBaseContact * contact,
									     const FrictionID & id,
									     const MemoryID & memory_id) :
  NTNFricLawCoulomb<Regularisation>(contact,id,memory_id),
  mu_s(0,1,0.,id+":mu_s",0.,"mu_s"),
  mu_k(0,1,0.,id+":mu_k",0.,"mu_k"),
  d_c(0,1,0.,id+":d_c",0.,"d_c") {
  AKANTU_DEBUG_IN();

  NTNFricLawCoulomb<Regularisation>::registerSynchronizedArray(this->mu_s);
  NTNFricLawCoulomb<Regularisation>::registerSynchronizedArray(this->mu_k);
  NTNFricLawCoulomb<Regularisation>::registerSynchronizedArray(this->d_c);

  this->registerParam("mu_s", this->mu_s, _pat_parsmod, "static friction coefficient");
  this->registerParam("mu_k", this->mu_k, _pat_parsmod, "kinetic friction coefficient");
  this->registerParam("d_c",  this->d_c,  _pat_parsmod, "slip weakening length");

  this->setParamAccessType("mu", _pat_readable);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawLinearSlipWeakening<Regularisation>::computeFrictionalStrength() {
  AKANTU_DEBUG_IN();

  computeFrictionCoefficient();
  NTNFricLawCoulomb<Regularisation>::computeFrictionalStrength();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawLinearSlipWeakening<Regularisation>::computeFrictionCoefficient() {
  AKANTU_DEBUG_IN();

  // get arrays
  const SynchronizedArray<bool> & stick = this->internalGetIsSticking();
  const SynchronizedArray<Real> & slip  = this->internalGetSlip();

  UInt nb_contact_nodes = this->contact->getNbContactNodes();
  for (UInt n=0; n<nb_contact_nodes; ++n) {
    if (stick(n)) {
      this->mu(n) = this->mu_s(n);
    }
    else {
      if (slip(n) >= this->d_c(n)) {
	this->mu(n) = this->mu_k(n);
      }
      else {
	// mu = mu_k + (1 - slip / Dc) * (mu_s - mu_k)
	this->mu(n) = this->mu_k(n)
		    + (1 - (slip(n) / this->d_c(n))) * (this->mu_s(n) - this->mu_k(n));
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawLinearSlipWeakening<Regularisation>::registerSynchronizedArray(SynchronizedArrayBase & array) {
  AKANTU_DEBUG_IN();

  this->mu_s.registerDependingArray(array);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawLinearSlipWeakening<Regularisation>::dumpRestart(const std::string & file_name) const {
  AKANTU_DEBUG_IN();

  this->mu_s.dumpRestartFile(file_name);
  this->mu_k.dumpRestartFile(file_name);
  this->d_c.dumpRestartFile(file_name);

  NTNFricLawCoulomb<Regularisation>::dumpRestart(file_name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawLinearSlipWeakening<Regularisation>::readRestart(const std::string & file_name) {
  AKANTU_DEBUG_IN();

  this->mu_s.readRestartFile(file_name);
  this->mu_k.readRestartFile(file_name);
  this->d_c.readRestartFile(file_name);

  NTNFricLawCoulomb<Regularisation>::readRestart(file_name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawLinearSlipWeakening<Regularisation>::printself(std::ostream & stream, int indent) const {
  AKANTU_DEBUG_IN();
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "NTNFricLawLinearSlipWeakening [" << std::endl;
  NTNFricLawCoulomb<Regularisation>::printself(stream, ++indent);
  stream << space << "]" << std::endl;


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawLinearSlipWeakening<Regularisation>::addDumpFieldToDumper(const std::string & dumper_name,
							     const std::string & field_id) {
  AKANTU_DEBUG_IN();

#ifdef AKANTU_USE_IOHELPER
  //  const SynchronizedArray<UInt> * nodal_filter = &(this->contact->getSlaves());

  if(field_id == "mu_s") {
    this->internalAddDumpFieldToDumper(dumper_name,
				       field_id,
				       new dumper::NodalField<Real>(this->mu_s.getArray()));
  }
  else if(field_id == "mu_k") {
    this->internalAddDumpFieldToDumper(dumper_name,
				       field_id,
				       new dumper::NodalField<Real>(this->mu_k.getArray()));
  }
  else if(field_id == "d_c") {
    this->internalAddDumpFieldToDumper(dumper_name,
				       field_id,
				       new dumper::NodalField<Real>(this->d_c.getArray()));
  }
  else {
    NTNFricLawCoulomb<Regularisation>::addDumpFieldToDumper(dumper_name, field_id);
  }

#endif

  AKANTU_DEBUG_OUT();
}

__END_AKANTU__