/**
 * @file   ntn_friclaw_linear_slip_weakening_no_healing.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Tue Dec 02 2014
 * @date last modification: Fri Jan 22 2016
 *
 * @brief  linear slip weakening
 *
 * @section LICENSE
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory
 * (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
#ifndef __AST_NTN_FRICLAW_LINEAR_SLIP_WEAKENING_NO_HEALING_HH__
#define __AST_NTN_FRICLAW_LINEAR_SLIP_WEAKENING_NO_HEALING_HH__

/* -------------------------------------------------------------------------- */
#include "ntn_friclaw_linear_slip_weakening.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template <class Regularisation = NTNFricRegNoRegularisation>
class NTNFricLawLinearSlipWeakeningNoHealing : public NTNFricLawLinearSlipWeakening<Regularisation> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  NTNFricLawLinearSlipWeakeningNoHealing(NTNBaseContact * contact,
					 const FrictionID & id = "linear_slip_weakening_no_healing",
					 const MemoryID & memory_id = 0);
  virtual ~NTNFricLawLinearSlipWeakeningNoHealing() {};
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;
  
protected:
  /// computes the friction coefficient as a function of slip
  virtual void computeFrictionCoefficient();

  /* ------------------------------------------------------------------------ */
  /* Dumpable                                                                 */
  /* ------------------------------------------------------------------------ */

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:

};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

/// standard output stream operator
template <class Regularisation>
inline std::ostream & operator <<(std::ostream & stream, 
				  const NTNFricLawLinearSlipWeakeningNoHealing<Regularisation> & _this)
{
  _this.printself(stream);
  return stream;
}

__END_AKANTU__

#include "ntn_friclaw_linear_slip_weakening_no_healing_tmpl.hh"

#endif /* __AST_NTN_FRICLAW_LINEAR_SLIP_WEAKENING_NO_HEALING_HH__ */
