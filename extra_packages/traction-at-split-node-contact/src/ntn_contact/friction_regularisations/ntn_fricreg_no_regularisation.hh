/**
 * @file   ntn_fricreg_no_regularisation.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 *
 * @brief  regularisation that does nothing
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */

#ifndef __AST_NTN_FRICREG_NO_REGULARISATION_HH__
#define __AST_NTN_FRICREG_NO_REGULARISATION_HH__

/* -------------------------------------------------------------------------- */
// simtools
#include "ntn_base_friction.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
class NTNFricRegNoRegularisation : public NTNBaseFriction {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  NTNFricRegNoRegularisation(NTNBaseContact * contact,
                             const FrictionID & id = "no_regularisation",
                             const MemoryID & memory_id = 0);
  virtual ~NTNFricRegNoRegularisation(){};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// set to steady state for no regularisation -> do nothing
  virtual void setToSteadyState(){};

  virtual void registerSynchronizedArray(SynchronizedArrayBase & array);
  virtual void dumpRestart(const std::string & file_name) const;
  virtual void readRestart(const std::string & file_name);

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

protected:
  virtual void computeFrictionalContactPressure();

  /// compute frictional strength according to friction law
  virtual void computeFrictionalStrength(){};

  /* ------------------------------------------------------------------------ */
  /* Dumpable                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  virtual void addDumpFieldToDumper(const std::string & dumper_name,
                                    const std::string & field_id);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
protected:
  /// get the is_in_contact array
  virtual const SynchronizedArray<bool> & internalGetIsInContact() {
    return this->contact->getIsInContact();
  };

  /// get the contact pressure (the norm: scalar value)
  virtual const SynchronizedArray<Real> & internalGetContactPressure();

  /// get the frictional strength array
  virtual SynchronizedArray<Real> & internalGetFrictionalStrength() {
    return this->frictional_strength;
  };

  /// get the is_sticking array
  virtual SynchronizedArray<bool> & internalGetIsSticking() {
    return this->is_sticking;
  }

  /// get the slip array
  virtual SynchronizedArray<Real> & internalGetSlip() { return this->slip; }

  /// get the slip array
  virtual SynchronizedArray<Real> & internalGetCumulativeSlip() {
    return this->cumulative_slip;
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  // contact pressure (absolut value) for computation of friction
  SynchronizedArray<Real> frictional_contact_pressure;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "ntn_fricreg_no_regularisation_inline_impl.cc"

/// standard output stream operator
inline std::ostream & operator<<(std::ostream & stream,
                                 const NTNFricRegNoRegularisation & _this) {
  _this.printself(stream);
  return stream;
}

} // namespace akantu

#endif /* __AST_NTN_FRICREG_NO_REGULARISATION_HH__ */
