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

__BEGIN_AKANTU__

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
  /// compute relative tangential field (complet array)
  /// relative to master nodes
  virtual void computeRelativeTangentialField(const Array<Real> & field,
					      Array<Real> & rel_tang_field) const;

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

};

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const MIIASYMContact & _this)
{
  _this.printself(stream);
  return stream;
}

__END_AKANTU__

#endif /* __AST_MIIASYM_CONTACT_HH__ */
