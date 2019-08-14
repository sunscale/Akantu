/**
 * @file   aka_static_memory.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  Memory management
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
#include <sstream>
#include <stdexcept>

/* -------------------------------------------------------------------------- */
#include "aka_static_memory.hh"

/* -------------------------------------------------------------------------- */
namespace akantu {

bool StaticMemory::is_instantiated = false;
StaticMemory * StaticMemory::single_static_memory = nullptr;
UInt StaticMemory::nb_reference = 0;

/* -------------------------------------------------------------------------- */
StaticMemory & StaticMemory::getStaticMemory() {
  if (!single_static_memory) {
    single_static_memory = new StaticMemory();
    is_instantiated = true;
  }

  nb_reference++;

  return *single_static_memory;
}

/* -------------------------------------------------------------------------- */
void StaticMemory::destroy() {
  nb_reference--;
  if (nb_reference == 0) {
    delete single_static_memory;
  }
}

/* -------------------------------------------------------------------------- */
StaticMemory::~StaticMemory() {
  AKANTU_DEBUG_IN();

  MemoryMap::iterator memory_it;
  for (memory_it = memories.begin(); memory_it != memories.end(); ++memory_it) {
    ArrayMap::iterator vector_it;
    for (vector_it = (memory_it->second).begin();
         vector_it != (memory_it->second).end(); ++vector_it) {
      delete vector_it->second;
    }
    (memory_it->second).clear();
  }
  memories.clear();
  is_instantiated = false;
  StaticMemory::single_static_memory = nullptr;
  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
void StaticMemory::sfree(const MemoryID & memory_id, const ID & name) {
  AKANTU_DEBUG_IN();

  try {
    auto & vectors = const_cast<ArrayMap &>(getMemory(memory_id));
    ArrayMap::iterator vector_it;
    vector_it = vectors.find(name);
    if (vector_it != vectors.end()) {
      AKANTU_DEBUG_INFO("Array " << name
                                 << " removed from the static memory number "
                                 << memory_id);
      delete vector_it->second;
      vectors.erase(vector_it);
      AKANTU_DEBUG_OUT();
      return;
    }
  } catch (debug::Exception & e) {
    AKANTU_EXCEPTION("The memory "
                     << memory_id << " does not exist (perhaps already freed) ["
                     << e.what() << "]");
    AKANTU_DEBUG_OUT();
    return;
  }

  AKANTU_DEBUG_INFO("The vector " << name
                                  << " does not exist (perhaps already freed)");
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void StaticMemory::printself(std::ostream & stream, int indent) const {
  std::string space = "";
  for (Int i = 0; i < indent; i++, space += AKANTU_INDENT)
    ;

  std::streamsize prec = stream.precision();
  stream.precision(2);

  stream << space << "StaticMemory [" << std::endl;
  UInt nb_memories = memories.size();
  stream << space << " + nb memories : " << nb_memories << std::endl;

  Real tot_size = 0;
  MemoryMap::const_iterator memory_it;
  for (memory_it = memories.begin(); memory_it != memories.end(); ++memory_it) {
    Real mem_size = 0;

    stream << space << AKANTU_INDENT << "Memory [" << std::endl;
    UInt mem_id = memory_it->first;
    stream << space << AKANTU_INDENT << " + memory id   : " << mem_id
           << std::endl;
    UInt nb_vectors = (memory_it->second).size();
    stream << space << AKANTU_INDENT << " + nb vectors  : " << nb_vectors
           << std::endl;
    stream.precision(prec);
    ArrayMap::const_iterator vector_it;
    for (vector_it = (memory_it->second).begin();
         vector_it != (memory_it->second).end(); ++vector_it) {
      (vector_it->second)->printself(stream, indent + 2);
      mem_size += (vector_it->second)->getMemorySize();
    }
    stream << space << AKANTU_INDENT
           << " + total size  : " << printMemorySize<char>(mem_size)
           << std::endl;
    stream << space << AKANTU_INDENT << "]" << std::endl;
    tot_size += mem_size;
  }
  stream << space << " + total size  : " << printMemorySize<char>(tot_size)
         << std::endl;
  stream << space << "]" << std::endl;

  stream.precision(prec);
}

/* -------------------------------------------------------------------------- */

} // namespace akantu
