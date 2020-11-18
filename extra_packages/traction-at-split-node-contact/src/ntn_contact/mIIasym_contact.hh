/**
 * @file   mIIasym_contact.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Tue Dec 02 2014
 * @date last modification: Fri Feb 23 2018
 *
 * @brief  contact for mode II anti-symmetric simulations
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

#ifndef AST_MIIASYM_CONTACT_HH_
#define AST_MIIASYM_CONTACT_HH_

/* -------------------------------------------------------------------------- */
// simtools
#include "ntrf_contact.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
class MIIASYMContact : public NTRFContact {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MIIASYMContact(SolidMechanicsModel & model, const ID & id = "contact",
                 const MemoryID & memory_id = 0);
  virtual ~MIIASYMContact() = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// update the impedance matrix
  virtual void updateImpedance();

  /// compute contact pressure -> do nothing because can only compute it in
  /// equilibrium
  virtual void computeContactPressure(){};

  /// compute relative normal field (only value that has to be multiplied with
  /// the normal)
  /// WARNING: this is only valid for the acceleration in equilibrium
  virtual void computeRelativeNormalField(const Array<Real> & field,
                                          Array<Real> & rel_normal_field) const;

  /// compute relative tangential field (complet array)
  /// relative to master nodes
  virtual void
  computeRelativeTangentialField(const Array<Real> & field,
                                 Array<Real> & rel_tang_field) const;

  /// compute contact pressure that is used over the entire time
  virtual void computeContactPressureInEquilibrium();

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;
};

/// standard output stream operator
inline std::ostream & operator<<(std::ostream & stream,
                                 const MIIASYMContact & _this) {
  _this.printself(stream);
  return stream;
}

} // namespace akantu

#endif /* AST_MIIASYM_CONTACT_HH_ */
