/**
 * @file   ntrf_friction.hh
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Mon Nov  5 10:21:11 2012
 *
 * @brief  friction for node to rigid flat interface
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

#ifndef __AST_NTRF_FRICTION_HH__
#define __AST_NTRF_FRICTION_HH__

/* -------------------------------------------------------------------------- */
// simtools
#include "ntrf_contact.hh"

__BEGIN_SIMTOOLS__

/* -------------------------------------------------------------------------- */
using namespace akantu;

/* -------------------------------------------------------------------------- */
class NTRFFriction : protected Memory, public Dumpable {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  NTRFFriction(NTRFContact & contact,
	       const FrictionID & id = "friction",
	       const MemoryID & memory_id = 0);
  virtual ~NTRFFriction() {};
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// compute friction traction
  virtual void computeFrictionTraction();
  
  /// compute stick traction (friction traction needed to stick the nodes)
  void computeStickTraction();

  /// apply the friction force
  void applyFrictionTraction();

  /// compute slip
  virtual void updateSlip();

  /// register syncronizedarrays for sync
  virtual void registerSyncronizedArray(SyncronizedArrayBase & array);
  
  /// dump restart file
  virtual void dumpRestart(const std::string & file_name) const;

  /// read restart file
  virtual void readRestart(const std::string & file_name);

  // set to steady state 
  virtual void setToSteadyState() {
    AKANTU_DEBUG_TO_IMPLEMENT();
  };

  /// get the number of sticking nodes (in parallel)
  virtual UInt getNbStickingNodes() const;

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

protected:
  /// compute frictional strength according to friction law
  virtual void computeFrictionalStrength() = 0;

  /* ------------------------------------------------------------------------ */
  /* Dumpable                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  virtual void addDumpFieldToDumper(const std::string & dumper_name,
				    const std::string & field_id);
  //  virtual void addDumpFieldVector(const std::string & field_id);
  // void dump() { this->contact.dump(); };

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(Contact, contact, const NTRFContact &)

  AKANTU_GET_MACRO(IsSticking,                 is_sticking, const SyncronizedArray<bool> &)
  AKANTU_GET_MACRO(FrictionalStrength, frictional_strength, const SyncronizedArray<Real> &)
  AKANTU_GET_MACRO(FrictionTraction,     friction_traction, const SyncronizedArray<Real> &)
  AKANTU_GET_MACRO(Slip,                              slip, const SyncronizedArray<Real> &)

protected:
  void setInternalArray(SyncronizedArray<Real> & array, Real value);
  void setInternalArray(SyncronizedArray<Real> & array, UInt node, Real value);
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  FrictionID id;

  NTRFContact & contact;
  
  // if node is sticking
  SyncronizedArray<bool> is_sticking;
  // frictional strength
  SyncronizedArray<Real> frictional_strength;
  // friction force
  SyncronizedArray<Real> friction_traction;
  // slip
  SyncronizedArray<Real> slip;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "ntrf_friction_inline_impl.cc"

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const NTRFFriction & _this)
{
  _this.printself(stream);
  return stream;
}

__END_SIMTOOLS__

#endif /* __AST_NTRF_FRICTION_HH__ */
