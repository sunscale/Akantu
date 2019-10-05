/**
 * @file   aka_memory.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  static memory wrapper
 *
 * @section LICENSE
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
#include <utility>

#include "aka_memory.hh"
#include "aka_static_memory.hh"

/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
Memory::Memory(ID id, MemoryID memory_id)
    : static_memory(StaticMemory::getStaticMemory()), id(std::move(id)),
      memory_id(memory_id) {}

/* -------------------------------------------------------------------------- */
Memory::~Memory() {
  if (StaticMemory::isInstantiated()) {
    std::list<ID>::iterator it;
    for (it = handeld_vectors_id.begin(); it != handeld_vectors_id.end();
         ++it) {
      AKANTU_DEBUG(dblAccessory, "Deleting the vector " << *it);
      static_memory.sfree(memory_id, *it);
    }
    static_memory.destroy();
  }

  handeld_vectors_id.clear();
}

/* -------------------------------------------------------------------------- */

} // namespace akantu
