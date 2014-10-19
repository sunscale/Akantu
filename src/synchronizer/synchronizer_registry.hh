/**
 * @file   synchronizer_registry.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date creation: Thu Jun 16 2011
 * @date last modification: Tue Nov 06 2012
 *
 * @brief  Registry of synchronizers
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */


#ifndef __AKANTU_SYNCHRONIZER_REGISTRY_HH__
#define __AKANTU_SYNCHRONIZER_REGISTRY_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "data_accessor.hh"
#include "synchronizer.hh"
/* -------------------------------------------------------------------------- */


__BEGIN_AKANTU__

class SynchronizerRegistry {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  SynchronizerRegistry(DataAccessor & data_accessor);
  virtual ~SynchronizerRegistry();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// synchronize operation
  void synchronize(SynchronizationTag tag);

  /// asynchronous synchronization
  void asynchronousSynchronize(SynchronizationTag tag);

  /// wait end of asynchronous synchronization
  void waitEndSynchronize(SynchronizationTag tag);

  /// register a new synchronization
  void registerSynchronizer(Synchronizer & synchronizer,SynchronizationTag tag);

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  // /// number of tags registered
  // UInt nb_synchronization_tags;

  typedef std::multimap<SynchronizationTag, Synchronizer *> Tag2Sync;
  /// list of registered synchronization
  Tag2Sync synchronizers;

  /// data accessor that will permit to do the pack/unpack things
  DataAccessor & data_accessor;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

// #include "synchronizer_registry_inline_impl.cc"

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const SynchronizerRegistry & _this)
{
  _this.printself(stream);
  return stream;
}


__END_AKANTU__



#endif /* __AKANTU_SYNCHRONIZER_REGISTRY_HH__ */
