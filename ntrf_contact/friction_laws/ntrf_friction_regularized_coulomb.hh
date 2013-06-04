/**
 * @file   ntrf_friction_regularized_coulomb.hh
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Thu May 23 15:43:40 2013
 *
 * @brief  regularization of contact pressure in coulomb friction
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
#ifndef __AST_NTRF_FRICTION_REGULARIZED_COULOMB_HH__
#define __AST_NTRF_FRICTION_REGULARIZED_COULOMB_HH__

/* -------------------------------------------------------------------------- */
// simtools
#include "ntrf_friction_coulomb.hh"

__BEGIN_SIMTOOLS__

/* -------------------------------------------------------------------------- */
using namespace akantu;

class NTRFFrictionRegularizedCoulomb : public NTRFFrictionCoulomb {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  NTRFFrictionRegularizedCoulomb(NTRFContact & contact,
				 const FrictionID & id = "friction_regularized_coulomb",
				 const MemoryID & memory_id = 0);
  virtual ~NTRFFrictionRegularizedCoulomb() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// register syncronizedarrays for sync
  virtual void registerSyncronizedArray(SyncronizedArrayBase & array);

  /// dump restart file
  virtual void dumpRestart(const std::string & file_name) const;

  /// read restart file
  virtual void readRestart(const std::string & file_name);

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

  virtual void setToSteadyState();

protected:
  /// compute the frictional contact pressure, what is used to compute friction = \mu * \fric_cont_pres
  virtual void computeFrictionalContactPressure();

  /* ------------------------------------------------------------------------ */
  /* Dumpable                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  virtual void addDumpField(const std::string & field_id);
  //  virtual void addDumpFieldVector(const std::string & field_id);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_SET_MACRO(RegularizationOn, regularization_on, bool);

  // set t_star to all nodes
  void setTStar(Real tstar);

  // set t_star only to node (global index)
  void setTStar(UInt node, Real tstar);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  // define if regularization is applied
  bool regularization_on;

  // characteristic time scale for regularisation of frict. contact pressure
  SyncronizedArray<Real> t_star;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "ntrf_friction_regularized_coulomb_inline_impl.cc"

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream,
				  const NTRFFrictionRegularizedCoulomb & _this)
{
  _this.printself(stream);
  return stream;
}

__END_SIMTOOLS__

#endif /* __AST_NTRF_FRICTION_REGULARIZED_COULOMB_HH__ */
