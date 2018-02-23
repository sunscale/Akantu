/**
 * @file   mIIasym_contact.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 *
 * @brief  contact for mode II anti-symmetric simulations
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */

#ifndef __AST_MIIASYM_CONTACT_HH__
#define __AST_MIIASYM_CONTACT_HH__

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

  MIIASYMContact(SolidMechanicsModel & model,
		 const ContactID & id = "contact",
		 const MemoryID & memory_id = 0);
  virtual ~MIIASYMContact() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// update the impedance matrix
  virtual void updateImpedance();

  /// compute contact pressure -> do nothing because can only compute it in equilibrium
  virtual void computeContactPressure() {};

  /// compute relative normal field (only value that has to be multiplied with the normal)
  /// WARNING: this is only valid for the acceleration in equilibrium
  virtual void computeRelativeNormalField(const Array<Real> & field,
					  Array<Real> & rel_normal_field) const;

  /// compute relative tangential field (complet array)
  /// relative to master nodes
  virtual void computeRelativeTangentialField(const Array<Real> & field,
					      Array<Real> & rel_tang_field) const;

  /// compute contact pressure that is used over the entire time
  virtual void computeContactPressureInEquilibrium();

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

};

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const MIIASYMContact & _this)
{
  _this.printself(stream);
  return stream;
}

} // namespace akantu

#endif /* __AST_MIIASYM_CONTACT_HH__ */
