/**
 * @file   ntn_friction.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 *
 * @brief  implementation of friction for node to node contact
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */

#ifndef __AST_NTN_FRICTION_HH__
#define __AST_NTN_FRICTION_HH__

/* -------------------------------------------------------------------------- */
// simtools
#include "ntn_base_friction.hh"
#include "ntn_friclaw_coulomb.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <template<class> class FrictionLaw = NTNFricLawCoulomb, 
	  class Regularisation = NTNFricRegNoRegularisation>
class NTNFriction : public FrictionLaw<Regularisation> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  NTNFriction(NTNBaseContact * contact,
	      const FrictionID & id = "friction",
	      const MemoryID & memory_id = 0);
  virtual ~NTNFriction() {};
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// apply the friction force
  virtual void applyFrictionTraction();

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

protected:

  /* ------------------------------------------------------------------------ */
  /* Dumpable                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  // virtual void addDumpFieldToDumper(const std::string & dumper_name,
  // 				    const std::string & field_id);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:

};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

/// standard output stream operator
template <template<class> class FrictionLaw, class Regularisation>
inline std::ostream & operator <<(std::ostream & stream, 
				  const NTNFriction<FrictionLaw,
				                    Regularisation> & _this)
{
  _this.printself(stream);
  return stream;
}

} // namespace akantu

#include "ntn_friction_tmpl.hh"

#endif /* __AST_NTN_FRICTION_HH__ */
