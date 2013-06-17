/**
 * @file   ntrf_friction_coulomb.hh
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Thu Mar 14 14:27:26 2013
 *
 * @brief  coulomb friction with \mu_s = \mu_k (constant)
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
#ifndef __AST_NTRF_FRICTION_COULOMB_HH__
#define __AST_NTRF_FRICTION_COULOMB_HH__

/* -------------------------------------------------------------------------- */
// simtools
#include "ntrf_friction.hh"

__BEGIN_SIMTOOLS__

/* -------------------------------------------------------------------------- */
using namespace akantu;

class NTRFFrictionCoulomb : public NTRFFriction {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  NTRFFrictionCoulomb(NTRFContact & contact,
		      const FrictionID & id = "friction_coulomb",
		      const MemoryID & memory_id = 0);
  virtual ~NTRFFrictionCoulomb() {};
  
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

  virtual void setToSteadyState() {};

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;
  
protected:
  /// compute frictional strength according to friction law
  virtual void computeFrictionalStrength();
  /// compute the frictional contact pressure, what is used to compute friction = \mu * \fric_cont_pres
  virtual void computeFrictionalContactPressure();

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
  // set friction coefficient to all nodes
  void setMu(Real mu);

  // set friction coefficient only to node (global index)
  void setMu(UInt node, Real mu);
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  // friction coefficient
  SynchronizedArray<Real> mu;

  // contact pressure (absolut value) for computation of friction
  SynchronizedArray<Real> frictional_contact_pressure;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "ntrf_friction_coulomb_inline_impl.cc"

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const NTRFFrictionCoulomb & _this)
{
  _this.printself(stream);
  return stream;
}

__END_SIMTOOLS__

#endif /* __AST_NTRF_FRICTION_COULOMB_HH__ */
