/**
 * @file   aka_static_memory_inline_impl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Jul 15 2010
 * @date last modification: Wed Feb 03 2016
 *
 * @brief  Implementation of inline functions of the class StaticMemory
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
#include "aka_static_memory.hh"

namespace akantu {
/* -------------------------------------------------------------------------- */
inline const ArrayMap &
StaticMemory::getMemory(const MemoryID & memory_id) const {
  AKANTU_DEBUG_IN();
  MemoryMap::const_iterator memory_it;
  memory_it = memories.find(memory_id);

  if (memory_it == memories.end()) {
    AKANTU_SILENT_EXCEPTION("StaticMemory as no memory with ID " << memory_id);
  }

  AKANTU_DEBUG_OUT();
  return memory_it->second;
}

/* -------------------------------------------------------------------------- */
inline const ArrayBase & StaticMemory::getArray(const MemoryID & memory_id,
                                                const ID & name) const {
  AKANTU_DEBUG_IN();

  const ArrayMap & vectors = getMemory(memory_id);

  ArrayMap::const_iterator vectors_it;
  vectors_it = vectors.find(name);
  if (vectors_it == vectors.end()) {
    AKANTU_SILENT_EXCEPTION("StaticMemory as no array named "
                            << name << " for the Memory " << memory_id);
  }

  AKANTU_DEBUG_OUT();
  return *(vectors_it->second);
}

/* -------------------------------------------------------------------------- */
} // namespace akantu
