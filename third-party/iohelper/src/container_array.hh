/**
 * @file   container_array.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date creation: Fri Oct 12 2012
 * @date last modification: Tue Feb 05 2013
 *
 * @brief  container array header
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
/* -------------------------------------------------------------------------- */
#include "iohelper_common.hh"
/* -------------------------------------------------------------------------- */

#ifndef IOHELPER_CONTAINER_ARRAY_HH_
#define IOHELPER_CONTAINER_ARRAY_HH_
/* -------------------------------------------------------------------------- */

namespace iohelper {
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
template <typename T>
class ContainerArray {
  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */

public:

  class iterator: public ::iohelper::iterator<T, iterator> {
  public:

    iterator(T * ptr, UInt dimension, UInt increment, const ElemType & el_type = MAX_ELEM_TYPE){
      this->ptr = ptr;
      this->increment = increment;
      this->dimension = dimension;
      this->el_type = el_type;
    };

    bool operator!=(const iterator & it) const override {
      return it.ptr != this->ptr;
    };

    iterator & operator++() override {
      ptr += increment;
      return *this;
    };

    IOHelperVector<T> operator*() override {
      return IOHelperVector<T>(ptr, increment);
    };

    ElemType element_type() override { return el_type; }

  private:

    T * ptr;
    UInt increment;
    UInt dimension;
    ElemType el_type;
  };


  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  ContainerArray(T * data, UInt dimension, UInt size, UInt stride=1,ElemType el_type = MAX_ELEM_TYPE){
    this->data = data;
    this->dimension = dimension;
    if (this->data == NULL) {
      this->_size = 0;
    } else {
      this->_size = size;
    }
    this->stride = stride;
    this->el_type = el_type;
  };

  virtual ~ContainerArray() = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */

  iterator begin(){
    return iterator(data,dimension,stride*dimension,el_type);
  };
  iterator end(){
    return iterator(data+_size*dimension*stride,dimension,stride*dimension);
  };

  UInt getDim(){return dimension;};
  UInt size(){return _size;};
  bool isHomogeneous(){return true;};

  DataType getDataType() { return ::iohelper::getDataType<T>(); }

  const ElemType & getElemType() { return el_type;}
  void setElemType(const ElemType & type) { el_type = type;}

public:

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  T * data;
  UInt dimension;
  UInt stride;
  UInt _size;
  ElemType el_type;
};

/* -------------------------------------------------------------------------- */


}

#endif /* IOHELPER_ITERATOR_ARRAY_HH_ */
