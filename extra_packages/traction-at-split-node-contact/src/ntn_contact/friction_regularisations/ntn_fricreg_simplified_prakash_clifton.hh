/**
 * @file   ntn_fricreg_simplified_prakash_clifton.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 *
 * @brief  regularisation that regularizes the frictional strength with one
 * parameter
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#ifndef __AST_NTN_FRICREG_SIMPLIFIED_PRAKASH_CLIFTON_HH__
#define __AST_NTN_FRICREG_SIMPLIFIED_PRAKASH_CLIFTON_HH__

/* -------------------------------------------------------------------------- */
// simtools
#include "ntn_fricreg_no_regularisation.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
class NTNFricRegSimplifiedPrakashClifton : public NTNFricRegNoRegularisation {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  NTNFricRegSimplifiedPrakashClifton(
      NTNBaseContact * contact,
      const FrictionID & id = "simplified_prakash_clifton",
      const MemoryID & memory_id = 0);
  virtual ~NTNFricRegSimplifiedPrakashClifton(){};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  virtual void registerSynchronizedArray(SynchronizedArrayBase & array);
  virtual void dumpRestart(const std::string & file_name) const;
  virtual void readRestart(const std::string & file_name);

  virtual void setToSteadyState();

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
public:
protected:
  /// get the frictional strength array
  virtual SynchronizedArray<Real> & internalGetFrictionalStrength() {
    return this->spc_internal;
  };

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  SynchronizedArray<Real> t_star;

  // to get the incremental frictional strength
  SynchronizedArray<Real> spc_internal;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "ntn_fricreg_simplified_prakash_clifton_inline_impl.cc"

/// standard output stream operator
inline std::ostream &
operator<<(std::ostream & stream,
           const NTNFricRegSimplifiedPrakashClifton & _this) {
  _this.printself(stream);
  return stream;
}

} // namespace akantu

#endif /* __AST_NTN_FRICREG_SIMPLIFIED_PRAKASH_CLIFTON_HH__ */
