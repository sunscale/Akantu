/**
 * @file   ntn_fricreg_simplified_prakash_clifton.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 *
 * @brief  implementation of simplified prakash clifton with one parameter
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
// simtools
#include "ntn_fricreg_simplified_prakash_clifton.hh"
#include "dumper_text.hh"
#include "dumper_nodal_field.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
NTNFricRegSimplifiedPrakashClifton::NTNFricRegSimplifiedPrakashClifton(NTNBaseContact * contact,
								       const FrictionID & id,
								       const MemoryID & memory_id) :
  NTNFricRegNoRegularisation(contact, id, memory_id),
  t_star(0,1,0.,id+":t_star",0.,"t_star"),
  spc_internal(0,1,0.,id+":spc_internal",0.,"spc_internal") {
  AKANTU_DEBUG_IN();
  
  NTNFricRegNoRegularisation::registerSynchronizedArray(this->t_star);
  NTNFricRegNoRegularisation::registerSynchronizedArray(this->spc_internal);

  this->registerParam("t_star", this->t_star, _pat_parsmod, "time scale of regularisation");
  this->registerParam("spc_internal", this->spc_internal, _pat_internal, "");

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFricRegSimplifiedPrakashClifton::computeFrictionalStrength() {
  AKANTU_DEBUG_IN();
  
  SolidMechanicsModel & model = this->contact->getModel();
  Real delta_t = model.getTimeStep();
  
  UInt nb_contact_nodes = this->contact->getNbContactNodes();
  for (UInt n=0; n<nb_contact_nodes; ++n) {
    Real alpha = delta_t / this->t_star(n);
    this->frictional_strength(n) += alpha * this->spc_internal(n);
    this->frictional_strength(n) /= 1 + alpha;
  }  

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFricRegSimplifiedPrakashClifton::setToSteadyState() {
  AKANTU_DEBUG_IN();

  /// fill the spc_internal array
  computeFrictionalStrength();

  /// set strength without regularisation
  UInt nb_contact_nodes = this->contact->getNbContactNodes();
  for (UInt n=0; n<nb_contact_nodes; ++n) {
    this->frictional_strength(n) = this->spc_internal(n);
  }  

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFricRegSimplifiedPrakashClifton::registerSynchronizedArray(SynchronizedArrayBase & array) {
  AKANTU_DEBUG_IN();
  
  this->t_star.registerDependingArray(array);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFricRegSimplifiedPrakashClifton::dumpRestart(const std::string & file_name) const {
  AKANTU_DEBUG_IN();
  
  this->t_star.dumpRestartFile(file_name);
  this->spc_internal.dumpRestartFile(file_name);

  NTNFricRegNoRegularisation::dumpRestart(file_name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFricRegSimplifiedPrakashClifton::readRestart(const std::string & file_name) {
  AKANTU_DEBUG_IN();
  
  this->t_star.readRestartFile(file_name);
  this->spc_internal.readRestartFile(file_name);

  NTNFricRegNoRegularisation::readRestart(file_name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFricRegSimplifiedPrakashClifton::printself(std::ostream & stream, 
						   int indent) const {
  AKANTU_DEBUG_IN();
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);
  
  stream << space << "NTNFricRegSimplifiedPrakashClifton [" << std::endl;
  NTNFricRegNoRegularisation::printself(stream, ++indent);
  stream << space << "]" << std::endl;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFricRegSimplifiedPrakashClifton::addDumpFieldToDumper(const std::string & dumper_name,
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

} // namespace akantu
