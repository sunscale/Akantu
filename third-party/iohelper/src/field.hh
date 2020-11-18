/**
 * @file   field.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date creation: Thu Mar 11 2010
 * @date last modification: Tue Feb 05 2013
 *
 * @brief  description of field
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef IOHELPER_FIELD_HH_
#define IOHELPER_FIELD_HH_
/* -------------------------------------------------------------------------- */
#include "field_interface.hh"
// #include "paraview_helper.hh"
// #include "dumper_lammps.hh"
/* -------------------------------------------------------------------------- */

namespace iohelper {

template <class Cont>
class Field : public FieldInterface {

  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */

public:
  using iterator = typename Cont::iterator;

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  Field(Cont & c, const std::string & name) : my_field(c), name(name){};
  ~Field() override = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  //! return true if the data is a constant size per element
  inline bool isHomogeneous() override;
  //! return the size per element (valid only if isHomogeneous is true)
  inline UInt getDim() override;
  //! return the number of stored items (elements, nodes, etc...)
  inline UInt size() override;
  //! return the description name of the container
  inline std::string getName() override;

  //! accept to be visited by a visitor
  void accept(Visitor & v) override;

  //! begin method
  iterator begin(){return my_field.begin();}
  //! end method
  iterator end(){return my_field.end();}


  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */

  inline DataType getDataType() override { return my_field.getDataType(); }

public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  Cont & my_field;

  std::string name;
};
/* -------------------------------------------------------------------------- */

template <class Cont>
bool Field<Cont>::isHomogeneous(){
  return my_field.isHomogeneous();
}

/* -------------------------------------------------------------------------- */


template <class Cont>
UInt Field<Cont>::getDim(){
  return my_field.getDim();
}

/* -------------------------------------------------------------------------- */


template <class Cont>
std::string Field<Cont>::getName(){
  return name;
}

/* -------------------------------------------------------------------------- */

template <class Cont>
UInt Field<Cont>::size(){
  return my_field.size();
}

/* -------------------------------------------------------------------------- */




}


#endif /* IOHELPER_FIELD_HH_ */
