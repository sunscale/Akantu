/**
 * @file   ntn_base_friction.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Tue Dec 02 2014
 * @date last modification: Fri Feb 23 2018
 *
 * @brief  base class for ntn and ntrf friction
 *
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#ifndef __AST_NTN_BASE_FRICTION_HH__
#define __AST_NTN_BASE_FRICTION_HH__

/* -------------------------------------------------------------------------- */
// akantu
#include "parsable.hh"

// simtools
#include "ntn_base_contact.hh"

/* -------------------------------------------------------------------------- */
namespace akantu {

/* -------------------------------------------------------------------------- */
template <>
inline void ParameterTyped<akantu::SynchronizedArray<Real>>::setAuto(
    const ParserParameter & in_param) {
  Parameter::setAuto(in_param);
  Real r = in_param;
  param.setAndChangeDefault(r);
}

/* -------------------------------------------------------------------------- */
template <>
template <>
inline void ParameterTyped<akantu::SynchronizedArray<Real>>::setTyped<Real>(
    const Real & value) {
  param.setAndChangeDefault(value);
}

/* -------------------------------------------------------------------------- */
class NTNBaseFriction : protected Memory, public Parsable, public Dumpable {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  NTNBaseFriction(NTNBaseContact & contact, const ID & id = "friction",
                  const MemoryID & memory_id = 0);
  virtual ~NTNBaseFriction() = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// compute friction traction
  virtual void computeFrictionTraction();

  /// compute stick traction (friction traction needed to stick the nodes)
  virtual void computeStickTraction();

  /// apply the friction force
  virtual void applyFrictionTraction();

  /// compute slip
  virtual void updateSlip();

  /// register Syncronizedarrays for sync
  virtual void registerSynchronizedArray(SynchronizedArrayBase & array);

  /// dump restart file
  virtual void dumpRestart(const std::string & file_name) const;

  /// read restart file
  virtual void readRestart(const std::string & file_name);

  /// set to steady state
  virtual void setToSteadyState() { AKANTU_TO_IMPLEMENT(); };

  /// get the number of sticking nodes (in parallel)
  /// a node that is not in contact does not count as sticking
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

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(Contact, contact, const NTNBaseContact &)

  AKANTU_GET_MACRO(IsSticking, is_sticking, const SynchronizedArray<bool> &)
  AKANTU_GET_MACRO(FrictionalStrength, frictional_strength,
                   const SynchronizedArray<Real> &)
  AKANTU_GET_MACRO(FrictionTraction, friction_traction,
                   const SynchronizedArray<Real> &)
  AKANTU_GET_MACRO(Slip, slip, const SynchronizedArray<Real> &)
  AKANTU_GET_MACRO(CumulativeSlip, cumulative_slip,
                   const SynchronizedArray<Real> &)
  AKANTU_GET_MACRO(SlipVelocity, slip_velocity, const SynchronizedArray<Real> &)

  /// set parameter of a given node
  /// (if you need to set to all: used the setMixed function of the Parsable).
  virtual void setParam(const std::string & name, UInt node, Real value);

  // replaced by the setMixed of the Parsable
  // virtual void setParam(const std::string & param, Real value) {
  //   AKANTU_ERROR("Friction does not know the following parameter: " <<
  //   param);
  // };

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  NTNBaseContact & contact;

  // if node is sticking
  SynchronizedArray<bool> is_sticking;
  // frictional strength
  SynchronizedArray<Real> frictional_strength;
  // friction force
  SynchronizedArray<Real> friction_traction;
  // slip
  SynchronizedArray<Real> slip;
  SynchronizedArray<Real> cumulative_slip;
  // slip velocity (tangential vector)
  SynchronizedArray<Real> slip_velocity;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "ntn_base_friction_inline_impl.hh"

/// standard output stream operator
inline std::ostream & operator<<(std::ostream & stream,
                                 const NTNBaseFriction & _this) {
  _this.printself(stream);
  return stream;
}

} // namespace akantu

#endif /* __AST_NTN_BASE_FRICTION_HH__ */
