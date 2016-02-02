/**
 * @file   ntn_fricreg_rubin_ampuero.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 *
 * @brief  implementation of no regularisation
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
// simtools
#include "ntn_fricreg_rubin_ampuero.hh"
#include "dumper_text.hh"
#include "dumper_nodal_field.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
NTNFricRegRubinAmpuero::NTNFricRegRubinAmpuero(NTNBaseContact * contact,
					       const FrictionID & id,
					       const MemoryID & memory_id) :
  NTNFricRegNoRegularisation(contact, id, memory_id),
  t_star(0,1,0.,id+":t_star",0.,"t_star") {
  AKANTU_DEBUG_IN();
  
  NTNFricRegNoRegularisation::registerSynchronizedArray(this->t_star);

  this->registerParam("t_star", this->t_star, _pat_parsmod, "time scale of regularization");

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
const SynchronizedArray<Real> & NTNFricRegRubinAmpuero::internalGetContactPressure() {
  AKANTU_DEBUG_IN();
  
  SolidMechanicsModel & model = this->contact->getModel();
  UInt dim = model.getSpatialDimension();
  Real delta_t = model.getTimeStep();
  
  // get contact arrays
  const SynchronizedArray<bool> & is_in_contact = this->internalGetIsInContact();
  const Array<Real> & pressure = this->contact->getContactPressure().getArray();
  Array<Real>::const_iterator< Vector<Real> > it = pressure.begin(dim);

  UInt nb_contact_nodes = this->contact->getNbContactNodes();
  for (UInt n=0; n<nb_contact_nodes; ++n) {
    // node pair is NOT in contact
    if (!is_in_contact(n))
      this->frictional_contact_pressure(n) = 0.;

    // if t_star is too small compute like Coulomb friction (without regularization)
    else if (Math::are_float_equal(this->t_star(n), 0.)) {
      const Vector<Real> & pres = it[n];
      this->frictional_contact_pressure(n) = pres.norm();
    }
    
    else {
      // compute frictional contact pressure
      // backward euler method: first order implicit numerical integration method
      // \reg_pres_n+1 = (\reg_pres_n + \delta_t / \t_star * \cur_pres)
      //               / (1 + \delta_t / \t_star)
      Real alpha = delta_t / this->t_star(n);
      const Vector<Real> & pres = it[n];
      this->frictional_contact_pressure(n) += alpha * pres.norm();
      this->frictional_contact_pressure(n) /= 1 + alpha;
    }
  }

  AKANTU_DEBUG_OUT();
  return this->frictional_contact_pressure;
}

/* -------------------------------------------------------------------------- */
void NTNFricRegRubinAmpuero::setToSteadyState() {
  AKANTU_DEBUG_IN();
  NTNFricRegNoRegularisation::computeFrictionalContactPressure();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFricRegRubinAmpuero::registerSynchronizedArray(SynchronizedArrayBase & array) {
  AKANTU_DEBUG_IN();
  
  this->t_star.registerDependingArray(array);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFricRegRubinAmpuero::dumpRestart(const std::string & file_name) const {
  AKANTU_DEBUG_IN();
  
  this->t_star.dumpRestartFile(file_name);

  NTNFricRegNoRegularisation::dumpRestart(file_name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFricRegRubinAmpuero::readRestart(const std::string & file_name) {
  AKANTU_DEBUG_IN();
  
  this->t_star.readRestartFile(file_name);

  NTNFricRegNoRegularisation::readRestart(file_name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFricRegRubinAmpuero::printself(std::ostream & stream, 
				       int indent) const {
  AKANTU_DEBUG_IN();
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);
  
  stream << space << "NTNFricRegRubinAmpuero [" << std::endl;
  NTNFricRegNoRegularisation::printself(stream, ++indent);
  stream << space << "]" << std::endl;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFricRegRubinAmpuero::addDumpFieldToDumper(const std::string & dumper_name,
						  const std::string & field_id) {
  AKANTU_DEBUG_IN();
  
#ifdef AKANTU_USE_IOHELPER
  //  const SynchronizedArray<UInt> * nodal_filter = &(this->contact->getSlaves());
  
  if (field_id == "t_star") {
    this->internalAddDumpFieldToDumper(dumper_name,
				       field_id,
				       new dumper::NodalField<Real>(this->t_star.getArray()));
  }
  else {
    NTNFricRegNoRegularisation::addDumpFieldToDumper(dumper_name, field_id);
  }
  
#endif

  AKANTU_DEBUG_OUT();
}

__END_AKANTU__