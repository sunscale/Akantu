/**
 * @file   mIIasym_contact.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 *
 * @brief  
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
// simtools
#include "mIIasym_contact.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
MIIASYMContact::MIIASYMContact(SolidMechanicsModel & model,
			       const ContactID & id,
			       const MemoryID & memory_id) :
  NTRFContact(model,id,memory_id)
{
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MIIASYMContact::updateImpedance() {
  AKANTU_DEBUG_IN();

  NTRFContact::updateImpedance();

  for (UInt i=0; i<this->impedance.getSize(); ++i) {
    this->impedance(i) *= 0.5;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/// WARNING: this is only valid for the acceleration in equilibrium
void MIIASYMContact::computeRelativeNormalField(const Array<Real> & field,
						Array<Real> & rel_normal_field) const {
  AKANTU_DEBUG_IN();

  NTRFContact::computeRelativeNormalField(field, rel_normal_field);

  for (Array<Real>::iterator<Real> it_rtfield  = rel_normal_field.begin();
       it_rtfield != rel_normal_field.end();
       ++it_rtfield) {

    // in the anti-symmetric case
    *it_rtfield *= 2.;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MIIASYMContact::computeRelativeTangentialField(const Array<Real> & field,
						    Array<Real> & rel_tang_field) const {
  AKANTU_DEBUG_IN();

  NTRFContact::computeRelativeTangentialField(field, rel_tang_field);

  UInt dim = this->model.getSpatialDimension();
  
  for (Array<Real>::iterator< Vector<Real> > it_rtfield  = rel_tang_field.begin(dim);
       it_rtfield != rel_tang_field.end(dim);
       ++it_rtfield) {

    // in the anti-symmetric case, the tangential fields become twice as large
    *it_rtfield *= 2.;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MIIASYMContact::computeContactPressureInEquilibrium() {
  AKANTU_DEBUG_IN();

  NTRFContact::computeContactPressure();

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void MIIASYMContact::printself(std::ostream & stream, int indent) const {
  AKANTU_DEBUG_IN();
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "MIIASYMContact [" << std::endl;
  NTRFContact::printself(stream,indent);
  stream << space << "]" << std::endl;

  AKANTU_DEBUG_OUT();
}

} // namespace akantu
