/**
 * @file   ntn_friction.hh
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Wed Feb 22 11:40:08 2012
 *
 * @brief  implementation of friction for node to node contact
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

#ifndef __AST_NTN_FRICTION_HH__
#define __AST_NTN_FRICTION_HH__

/* -------------------------------------------------------------------------- */
// simtools
#include "ntn_base_friction.hh"
#include "ntn_contact.hh"

__BEGIN_SIMTOOLS__

/* -------------------------------------------------------------------------- */
using namespace akantu;

/* -------------------------------------------------------------------------- */
class NTNFriction : public NTNBaseFriction {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  NTNFriction(NTNContact & contact,
	      const FrictionID & id = "friction",
	      const MemoryID & memory_id = 0);
  virtual ~NTNFriction() {};
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// apply the friction force
  virtual void applyFrictionTraction();

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

protected:
  /// compute frictional strength according to friction law
  virtual void computeFrictionalStrength() = 0;

  /* ------------------------------------------------------------------------ */
  /* Dumpable                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  // virtual void addDumpFieldToDumper(const std::string & dumper_name,
  // 				    const std::string & field_id);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:

};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const NTNFriction & _this)
{
  _this.printself(stream);
  return stream;
}

__END_SIMTOOLS__

#include "ntn_friction_tmpl.hh"

#endif /* __AST_NTN_FRICTION_HH__ */
