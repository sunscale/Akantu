/**
 * @file   variable.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Tue Jun 04 2013
 * @date last modification: Tue Jun 04 2013
 *
 * @brief  for dump of global variables
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
#ifndef IOHELPER_VARIABLE_HH_
#define IOHELPER_VARIABLE_HH_
/* -------------------------------------------------------------------------- */
#include "variable_interface.hh"
/* -------------------------------------------------------------------------- */

namespace iohelper {

template <class Cont>
class Variable : public VariableInterface {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  Variable(Cont & c, const std::string & name) : my_variable(c), name(name){};
  ~Variable() override = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  //! accept to be visited by a visitor
  inline void accept(Visitor & v) const override;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  inline DataType getDataType() const override {
    return my_variable.getDataType();
  }

  //! return the dim
  inline UInt getDim() const override { return my_variable.getDim(); }

  //! return the description name of the container
  inline std::string getName() const override { return name; }

  inline const Cont & operator*() const { return my_variable; } 

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  Cont & my_variable;

  std::string name;
};

/* -------------------------------------------------------------------------- */

}

#endif /* IOHELPER_VARIABLE_HH_ */
