/**
 * @file   ntn_friction_linear_slip_weakening.hh
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Tue Nov 20 14:56:20 2012
 *
 * @brief  linear slip weakening friction
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
#ifndef __AST_FRICTION_LINEAR_SLIP_WEAKENING_HH__
#define __AST_NTN_FRICTION_LINEAR_SLIP_WEAKENING_HH__

/* -------------------------------------------------------------------------- */
// simtools
#include "ntn_contact.hh"
#include "ntn_friction_coulomb.hh"

__BEGIN_SIMTOOLS__

/* -------------------------------------------------------------------------- */
using namespace akantu;

/* -------------------------------------------------------------------------- */
class NTNFrictionLinearSlipWeakening : public NTNFrictionCoulomb {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  NTNFrictionLinearSlipWeakening(NTNContact & contact,
				 const FrictionID & id = 
				 "friction_linear_slip_weakening",
				 const MemoryID & memory_id = 0);
  virtual ~NTNFrictionLinearSlipWeakening() {};
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// register syncronizedarrays for sync
  virtual void registerSyncronizedArray(SyncronizedArrayBase & array);
  
  /// dump restart file
  virtual void dumpRestart(std::string file_name) const;

  /// read restart file
  virtual void readRestart(std::string file_name);

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

protected:
  /// compute frictional strength according to friction law
  virtual void computeFrictionalStrength();

  // computes the friction coefficient as a function of slip
  virtual void computeFrictionCoefficient();

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  // set static friction coefficient to all nodes
  void setMuS(Real mu);

  // set static friction coefficient only to node (global index)
  void setMuS(UInt node, Real mu);

  // set kinetic friction coefficient to all nodes
  void setMuK(Real mu);

  // set kinetic friction coefficient only to node (global index)
  void setMuK(UInt node, Real mu);

  // set weakening length to all nodes
  void setWeakeningLength(Real length);

  // set weakening length only to node (global index)
  void setWeakeningLength(UInt node, Real length);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  // Dc the length over which slip weakening happens
  SyncronizedArray<Real> weakening_length;

  // static coefficient of friction
  SyncronizedArray<Real> mu_s;

  // kinetic coefficient of friction
  SyncronizedArray<Real> mu_k;

  // internal variable = absolut value of tangential gap when last sticked
  SyncronizedArray<Real> slip;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "ntn_friction_linear_slip_weakening_inline_impl.cc"

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const NTNFrictionLinearSlipWeakening & _this)
{
  _this.printself(stream);
  return stream;
}

__END_SIMTOOLS__

#endif /* __AST_NTN_FRICTION_LINEAR_SLIP_WEAKENING_HH__ */
