/**
 * @file   ntn_friclaw_linear_slip_weakening.hh
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Mon Sep  2 15:28:25 2013
 *
 * @brief  linear slip weakening
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
#ifndef __AST_NTN_FRICLAW_LINEAR_SLIP_WEAKENING_HH__
#define __AST_NTN_FRICLAW_LINEAR_SLIP_WEAKENING_HH__

/* -------------------------------------------------------------------------- */
// simtools
#include "ntn_friclaw_coulomb.hh"

__BEGIN_SIMTOOLS__

/* -------------------------------------------------------------------------- */
using namespace akantu;

template <class Regularisation = NTNFricRegNoRegularisation>
class NTNFricLawLinearSlipWeakening : public NTNFricLawCoulomb<Regularisation> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  NTNFricLawLinearSlipWeakening(NTNBaseContact * contact,
				const FrictionID & id = "linear_slip_weakening",
				const MemoryID & memory_id = 0);
  virtual ~NTNFricLawLinearSlipWeakening() {};
  
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
  /// computes the friction coefficient as a function of slip
  virtual void computeFrictionCoefficient();

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
  virtual void setParam(const std::string & param, UInt node, Real value);
  virtual void setParam(const std::string & param, Real value);
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  // Dc the length over which slip weakening happens
  SynchronizedArray<Real> d_c;

  // static coefficient of friction
  SynchronizedArray<Real> mu_s;

  // kinetic coefficient of friction
  SynchronizedArray<Real> mu_k;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

/// standard output stream operator
template <class Regularisation>
inline std::ostream & operator <<(std::ostream & stream, 
				  const NTNFricLawLinearSlipWeakening<Regularisation> & _this)
{
  _this.printself(stream);
  return stream;
}

__END_SIMTOOLS__

#include "ntn_friclaw_linear_slip_weakening_tmpl.hh"

#endif /* __AST_NTN_FRICLAW_LINEAR_SLIP_WEAKENING_HH__ */
