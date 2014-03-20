/**
 * @file   ntn_friclaw_linear_cohesive_tmpl.hh
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Thu Mar 20 15:23:39 2014
 *
 * @brief  implementation of linear cohesive law
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
//#include "dumper_text.hh"

__BEGIN_SIMTOOLS__

/* -------------------------------------------------------------------------- */
template <class Regularisation>
NTNFricLawLinearCohesive<Regularisation>::NTNFricLawLinearCohesive(NTNBaseContact * contact,
								   const FrictionID & id,
								   const MemoryID & memory_id) :
  Regularisation(contact,id,memory_id),
  G_c(0,1,0.,id+":G_c",0.,"G_c"),
  sigma_c(0,1,0.,id+":sigma_c",0.,"sigma_c") {
  AKANTU_DEBUG_IN();

  Regularisation::registerSynchronizedArray(this->G_c);
  Regularisation::registerSynchronizedArray(this->sigma_c);
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawLinearCohesive<Regularisation>::computeFrictionalStrength() {
  AKANTU_DEBUG_IN();

  SolidMechanicsModel & model = this->contact->getModel();
  UInt dim = model.getSpatialDimension();

  // get arrays
  const SynchronizedArray<bool> & is_in_contact = this->internalGetIsInContact();
  const SynchronizedArray<Real> & slip  = this->internalGetSlip();

  UInt nb_contact_nodes = this->contact->getNbContactNodes();
  for (UInt n=0; n<nb_contact_nodes; ++n) {
    // node pair is NOT in contact
    if (!is_in_contact(n))
      strength(n) = 0.;

    // node pair is in contact
    else {
      Real slope = (this->sigma_c(n) * this->sigma_c(n)) / (2*this->G_c(n));
      strength(n) = std::max(this->sigma_c(n) - slope * slip(n), 0.); // no negative strength
      // if we want to keep strength loss after restick, 
      // we need to do min between new strength and previous strength
    }
  }

  Regularisation::computeFrictionalStrength();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawLinearCohesive<Regularisation>::registerSynchronizedArray(SynchronizedArrayBase & array) {
  AKANTU_DEBUG_IN();
  
  this->G_c.registerDependingArray(array);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawLinearCohesive<Regularisation>::dumpRestart(const std::string & file_name) const {
  AKANTU_DEBUG_IN();
  
  this->G_c.dumpRestartFile(file_name);
  this->sigma_c.dumpRestartFile(file_name);

  Regularisation::dumpRestart(file_name);
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawLinearCohesive<Regularisation>::readRestart(const std::string & file_name) {
  AKANTU_DEBUG_IN();
  
  this->G_c.readRestartFile(file_name);
  this->sigma_c.readRestartFile(file_name);

  Regularisation::readRestart(file_name);
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawLinearCohesive<Regularisation>::setParam(const std::string & param, 
							Real value) {
  AKANTU_DEBUG_IN();

  if (param == "G_c") {
    this->setInternalArray(this->G_c, value);
  }
  else if (param == "sigma_c") {
    this->setInternalArray(this->sigma_c, value);
  }
  else {
    Regularisation::setParam(param, value);
  }
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawLinearCohesive<Regularisation>::setParam(const std::string & param, 
							UInt node, Real value) {
  AKANTU_DEBUG_IN();

  if (param == "G_c") {
    this->setInternalArray(this->G_c, node, value);
  }
  else if (param == "sigma_c") {
    this->setInternalArray(this->sigma_c, node, value);
  }
  else {
    Regularisation::setParam(param, value);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawLinearCohesive<Regularisation>::printself(std::ostream & stream,
							 int indent) const {
  AKANTU_DEBUG_IN();
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);
  
  stream << space << "NTNFricLawLinearCohesive [" << std::endl;
  stream << space << this->G_c << std::endl;
  stream << space << this->sigma_c << std::endl;
  Regularisation::printself(stream, ++indent);
  stream << space << "]" << std::endl;
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawLinearCohesive<Regularisation>::addDumpFieldToDumper(const std::string & dumper_name,
							     const std::string & field_id) {
  AKANTU_DEBUG_IN();
  
#ifdef AKANTU_USE_IOHELPER
  //  const SynchronizedArray<UInt> * nodal_filter = &(this->contact->getSlaves());
  
  if(field_id == "G_c") {
    this->internalAddDumpFieldToDumper(dumper_name,
				       field_id,
				       new DumperIOHelper::NodalField<Real>(this->G_c.getArray()));
  }
  else if(field_id == "sigma_c") {
    this->internalAddDumpFieldToDumper(dumper_name,
				       field_id,
				       new DumperIOHelper::NodalField<Real>(this->sigma_c.getArray()));
  }
  else {
    Regularisation::addDumpFieldToDumper(dumper_name, field_id);
  }
  
#endif

  AKANTU_DEBUG_OUT();
}

__END_SIMTOOLS__
