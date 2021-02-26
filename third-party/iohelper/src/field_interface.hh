/**
 * @file   field_interface.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date creation: Wed Dec 14 2011
 * @date last modification: Tue Feb 05 2013
 *
 * @brief  header for the field interface
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

#ifndef IOHELPER_CONTAINER_INTERFACE_HH_
#define IOHELPER_CONTAINER_INTERFACE_HH_
/* -------------------------------------------------------------------------- */
#include "visitor.hh"
/* -------------------------------------------------------------------------- */


namespace iohelper {

class FieldInterface {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  FieldInterface() = default;
  virtual ~FieldInterface() = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  virtual void accept(Visitor & v) = 0;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  //! return true if the data is a constant size per element
  virtual bool isHomogeneous()=0;
  //! return the size per element (valid only if isHomogeneous is true)
  virtual UInt getDim()=0;
  //! return the description name of the container
  virtual std::string getName()=0;
  //! return the number of stored items (elements, nodes, etc...)
  virtual UInt size()=0;

  virtual DataType getDataType()=0;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

};


/* -------------------------------------------------------------------------- */

}

#endif /* IOHELPER_CONTAINER_INTERFACE_HH_ */
