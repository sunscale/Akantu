/**
 * @file   aka_static_memory.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Sun Dec 03 2017
 *
 * @brief  Memory management
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
 * @section DESCRIPTION
 *
 * The class handling the memory, allocation/reallocation/desallocation
 * The objects can register their array and ask for allocation or realocation
 *
 */

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_STATIC_MEMORY_HH__
#define __AKANTU_STATIC_MEMORY_HH__

/* -------------------------------------------------------------------------- */
#include "aka_array.hh"
#include "aka_common.hh"

/* -------------------------------------------------------------------------- */
#include <map>

/* -------------------------------------------------------------------------- */
namespace akantu {

using ArrayMap = std::map<ID, ArrayBase *>;
using MemoryMap = std::map<MemoryID, ArrayMap>;

/**
 * @class StaticMemory
 * @brief Class for memory management common to all objects (this class as to
 * be accessed by an interface class memory)
 */
class StaticMemory {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
private:
  /// Default constructor
  StaticMemory() = default;

public:
  virtual ~StaticMemory();

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// Get the global instance of the StaticMemory
  static StaticMemory & getStaticMemory();

  static bool isInstantiated() { return is_instantiated; };

  /// remove a reference on the static memory
  void destroy();

  /// access to an Array
  inline const ArrayBase & getArray(const MemoryID & memory_id,
                                    const ID & name) const;

  /// get all vectors of a memory
  inline const ArrayMap & getMemory(const MemoryID & memory_id) const;

  /* ------------------------------------------------------------------------ */
  /* Class Methods                                                            */
  /* ------------------------------------------------------------------------ */
public:
  /**
   * Allocation of an array of type
   *
   * @param memory_id the id of the memory accessing to the static memory
   * @param name name of the array (for example connectivity)
   * @param size number of size (for example number of nodes)
   * @param nb_component number of component (for example spatial dimension)
   *
   * @return pointer an array of memory actual size: size * nb_component * sizeof(T)
   */
  template <typename T>
  Array<T> & smalloc(const MemoryID & memory_id, const ID & name, UInt size,
                     UInt nb_component);

  template <typename T>
  Array<T> & smalloc(const MemoryID & memory_id, const ID & name, UInt size,
                     UInt nb_component, const T & init_value);
  /**
   * free the memory associated to the array name
   *
   * @param memory_id the id of the memory accessing to the static memory
   * @param name the name of an existing array
   */
  void sfree(const MemoryID & memory_id, const ID & name);

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// is the static memory instantiated
  static bool is_instantiated;

  /// unique instance of the StaticMemory
  static StaticMemory * single_static_memory;

  /// map of all allocated arrays, indexed by their names
  MemoryMap memories;

  /// number of references on the static memory
  static UInt nb_reference;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

/// standard output stream operator
inline std::ostream & operator<<(std::ostream & stream,
                                 const StaticMemory & _this) {
  _this.printself(stream);
  return stream;
}

} // namespace akantu

#include "aka_static_memory_inline_impl.hh"
#include "aka_static_memory_tmpl.hh"

#endif /* __AKANTU_STATIC_MEMORY_HH__ */
