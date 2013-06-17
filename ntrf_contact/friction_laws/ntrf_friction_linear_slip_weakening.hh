/**
 * @file   ntrf_friction_linear_slip_weakening.hh
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Mon May 13 17:05:07 2013
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
#ifndef __AST_NTRF_FRICTION_LINEAR_SLIP_WEAKENING_HH__
#define __AST_NTRF_FRICTION_LINEAR_SLIP_WEAKENING_HH__

/* -------------------------------------------------------------------------- */
// simtools
#include "ntrf_friction_regularized_coulomb.hh"

__BEGIN_SIMTOOLS__

/* -------------------------------------------------------------------------- */
using namespace akantu;

/* -------------------------------------------------------------------------- */
class NTRFFrictionLinearSlipWeakening : public NTRFFrictionRegularizedCoulomb {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  NTRFFrictionLinearSlipWeakening(NTRFContact & contact,
				  const FrictionID & id = "friction_linear_slip_weakening",
				  const MemoryID & memory_id = 0);
  virtual ~NTRFFrictionLinearSlipWeakening() {};
  
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
  /// computes the friction coefficient as a function of slip
  virtual void computeFrictionCoefficient();
  /// compute frictional strength according to friction law
  virtual void computeFrictionalStrength();

  /* ------------------------------------------------------------------------ */
  /* Dumpable                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  virtual void addDumpFieldToDumper(const std::string & dumper_name,
				    const std::string & field_id);
  //  virtual void addDumpFieldVector(const std::string & field_id);
  
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
private:
  // Dc the length over which slip weakening happens
  SynchronizedArray<Real> weakening_length;

  // static coefficient of friction
  SynchronizedArray<Real> mu_s;

  // kinetic coefficient of friction
  SynchronizedArray<Real> mu_k;
 
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "ntrf_friction_linear_slip_weakening_inline_impl.cc"

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const NTRFFrictionLinearSlipWeakening & _this)
{
  _this.printself(stream);
  return stream;
}




__END_SIMTOOLS__

#endif /* __AST_NTRF_FRICTION_LINEAR_SLIP_WEAKENING_HH__ */
