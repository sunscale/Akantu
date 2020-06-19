/**
 * @file   aka_memory_inline_impl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Jul 15 2010
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  Implementation of the inline functions of the class Memory
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
#include "aka_memory.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {
/* -------------------------------------------------------------------------- */
template <class T>
inline Array<T> & Memory::alloc(const ID & name, UInt size, UInt nb_component) {
  handeld_vectors_id.push_back(name);
  return static_memory.smalloc<T>(memory_id, name, size, nb_component);
}

/* -------------------------------------------------------------------------- */
template <class T>
inline Array<T> & Memory::alloc(const ID & name, UInt size, UInt nb_component,
                                const T & init_value) {
  handeld_vectors_id.push_back(name);
  return static_memory.smalloc<T>(memory_id, name, size, nb_component,
                                  init_value);
}

/* -------------------------------------------------------------------------- */
inline void Memory::dealloc(const ID & name) {
  AKANTU_DEBUG(dblAccessory, "Deleting the vector " << name);
  static_memory.sfree(memory_id, name);
  handeld_vectors_id.remove(name);
}

/* -------------------------------------------------------------------------- */
template <class T> inline Array<T> & Memory::getArray(const ID & name) {
  return static_cast<Array<T> &>(
      const_cast<ArrayBase &>(static_memory.getArray(memory_id, name)));
}

/* -------------------------------------------------------------------------- */
template <class T>
inline const Array<T> & Memory::getArray(const ID & name) const {
  return static_cast<Array<T> &>(
      const_cast<ArrayBase &>(static_memory.getArray(memory_id, name)));
}

} // namespace akantu
