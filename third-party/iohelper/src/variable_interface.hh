/**
 * @file   variable_interface.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Tue Jun 04 2013
 * @date last modification: Tue Jun 04 2013
 *
 * @brief  variable interface
 *
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * IOHelper is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * IOHelper is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with IOHelper. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#ifndef IOHELPER_VARIABLE_INTERFACE_HH_
#define IOHELPER_VARIABLE_INTERFACE_HH_

/* -------------------------------------------------------------------------- */
#include "visitor.hh"

/* -------------------------------------------------------------------------- */
namespace iohelper {

class VariableInterface {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  VariableInterface() = default;
  virtual ~VariableInterface() = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  virtual void accept(Visitor & v) const = 0;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  //! return the size per element (valid only if isHomogeneous is true)
  virtual UInt getDim() const = 0;
  //! return the description name of the container
  virtual std::string getName() const = 0;

  virtual DataType getDataType() const = 0;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

};

/* -------------------------------------------------------------------------- */

}

#endif /* IOHELPER_VARIABLE_INTERFACE_HH_ */
