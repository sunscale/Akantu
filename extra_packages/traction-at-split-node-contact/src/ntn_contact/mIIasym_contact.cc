/**
 * @file   mIIasym_contact.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Tue Dec 02 2014
 * @date last modification: Fri Feb 23 2018
 *
 * @brief
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
#include "mIIasym_contact.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
MIIASYMContact::MIIASYMContact(SolidMechanicsModel & model, const ID & id,
                               const MemoryID & memory_id)
    : NTRFContact(model, id, memory_id) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MIIASYMContact::updateImpedance() {
  AKANTU_DEBUG_IN();

  NTRFContact::updateImpedance();

  for (UInt i = 0; i < this->impedance.size(); ++i) {
    this->impedance(i) *= 0.5;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/// WARNING: this is only valid for the acceleration in equilibrium
void MIIASYMContact::computeRelativeNormalField(
    const Array<Real> & field, Array<Real> & rel_normal_field) const {
  AKANTU_DEBUG_IN();

  NTRFContact::computeRelativeNormalField(field, rel_normal_field);

  for (auto it_rtfield = rel_normal_field.begin();
       it_rtfield != rel_normal_field.end(); ++it_rtfield) {

    // in the anti-symmetric case
    *it_rtfield *= 2.;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MIIASYMContact::computeRelativeTangentialField(
    const Array<Real> & field, Array<Real> & rel_tang_field) const {
  AKANTU_DEBUG_IN();

  NTRFContact::computeRelativeTangentialField(field, rel_tang_field);

  UInt dim = this->model.getSpatialDimension();

  for (Array<Real>::iterator<Vector<Real>> it_rtfield =
           rel_tang_field.begin(dim);
       it_rtfield != rel_tang_field.end(dim); ++it_rtfield) {

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
  for (Int i = 0; i < indent; i++, space += AKANTU_INDENT)
    ;

  stream << space << "MIIASYMContact [" << std::endl;
  NTRFContact::printself(stream, indent);
  stream << space << "]" << std::endl;

  AKANTU_DEBUG_OUT();
}

} // namespace akantu
