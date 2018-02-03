/**
 * @file   ntn_friclaw_coulomb.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 *
 * @brief  coulomb friction with \mu_s = \mu_k (constant)
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#ifndef __AST_NTN_FRICLAW_COULOMB_HH__
#define __AST_NTN_FRICLAW_COULOMB_HH__

/* -------------------------------------------------------------------------- */
// simtools
#include "ntn_fricreg_no_regularisation.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <class Regularisation = NTNFricRegNoRegularisation>
class NTNFricLawCoulomb : public Regularisation {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  NTNFricLawCoulomb(NTNBaseContact * contact,
		    const FrictionID & id = "coulomb",
		    const MemoryID & memory_id = 0);
  virtual ~NTNFricLawCoulomb() {};
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// register synchronizedarrays for sync
  virtual void registerSynchronizedArray(SynchronizedArrayBase & array);

  /// dump restart file
  virtual void dumpRestart(const std::string & file_name) const;

  /// read restart file
  virtual void readRestart(const std::string & file_name);

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;
  
protected:
  /// compute frictional strength according to friction law
  virtual void computeFrictionalStrength();

  /* ------------------------------------------------------------------------ */
  /* Dumpable                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  virtual void addDumpFieldToDumper(const std::string & dumper_name,
				    const std::string & field_id);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  // friction coefficient
  SynchronizedArray<Real> mu;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

/// standard output stream operator
template <class Regularisation>
inline std::ostream & operator <<(std::ostream & stream, 
				  const NTNFricLawCoulomb<Regularisation> & _this)
{
  _this.printself(stream);
  return stream;
}

} // namespace akantu

#include "ntn_friclaw_coulomb_tmpl.hh"

#endif /* __AST_NTN_FRICLAW_COULOMB_HH__ */
