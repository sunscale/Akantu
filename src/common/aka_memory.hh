/**
 * @file   aka_memory.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  wrapper for the static_memory, all object which wants
 * to access the static_memory as to inherit from the class memory
 *
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_array.hh"
#include "aka_static_memory.hh"
/* -------------------------------------------------------------------------- */
#include <list>
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
#ifndef __AKANTU_MEMORY_HH__
#define __AKANTU_MEMORY_HH__

namespace akantu {

class Memory {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
protected:
  Memory(ID id, MemoryID memory_id = 0);

  virtual ~Memory();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
protected:
  /// malloc
  template <class T>
  inline Array<T> & alloc(const ID & name, UInt size, UInt nb_component);

  /// malloc
  template <class T>
  inline Array<T> & alloc(const ID & name, UInt size, UInt nb_component,
                          const T & init_value);

  /* ------------------------------------------------------------------------ */
  /// free an array
  inline void dealloc(const ID & name);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
protected:
  template <typename T> inline Array<T> & getArray(const ID & name);

  template <typename T> inline const Array<T> & getArray(const ID & name) const;

public:
  AKANTU_GET_MACRO(MemoryID, memory_id, const MemoryID &);

  AKANTU_GET_MACRO(ID, id, const ID &);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// the static memory instance
  StaticMemory & static_memory;

  /// list of allocated vectors id
  std::list<ID> handeld_vectors_id;

protected:
  ID id;

  /// the id registred in the static memory
  MemoryID memory_id;
};

/* -------------------------------------------------------------------------- */
/* Inline functions                                                           */
/* -------------------------------------------------------------------------- */
} // namespace akantu

#include "aka_memory_inline_impl.hh"

#endif /* __AKANTU_MEMORY_HH__ */
