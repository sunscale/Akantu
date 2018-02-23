/**
 * @file   ntn_friclaw_linear_slip_weakening_no_healing.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 *
 * @brief  linear slip weakening
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#ifndef __AST_NTN_FRICLAW_LINEAR_SLIP_WEAKENING_NO_HEALING_HH__
#define __AST_NTN_FRICLAW_LINEAR_SLIP_WEAKENING_NO_HEALING_HH__

/* -------------------------------------------------------------------------- */
#include "ntn_friclaw_linear_slip_weakening.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <class Regularisation = NTNFricRegNoRegularisation>
class NTNFricLawLinearSlipWeakeningNoHealing
    : public NTNFricLawLinearSlipWeakening<Regularisation> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  NTNFricLawLinearSlipWeakeningNoHealing(
      NTNBaseContact * contact,
      const FrictionID & id = "linear_slip_weakening_no_healing",
      const MemoryID & memory_id = 0);
  virtual ~NTNFricLawLinearSlipWeakeningNoHealing(){};

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
inline std::ostream & operator<<(
    std::ostream & stream,
    const NTNFricLawLinearSlipWeakeningNoHealing<Regularisation> & _this) {
  _this.printself(stream);
  return stream;
}

} // namespace akantu

#include "ntn_friclaw_linear_slip_weakening_no_healing_tmpl.hh"

#endif /* __AST_NTN_FRICLAW_LINEAR_SLIP_WEAKENING_NO_HEALING_HH__ */
