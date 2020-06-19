/**
 * @file   synchronizer_registry.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Jun 16 2011
 * @date last modification: Wed Nov 15 2017
 *
 * @brief  Registry of synchronizers
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
#include "synchronizer_registry.hh"
#include "synchronizer.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
SynchronizerRegistry::SynchronizerRegistry() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
SynchronizerRegistry::~SynchronizerRegistry() {
  AKANTU_DEBUG_IN();

  synchronizers.clear();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SynchronizerRegistry::registerDataAccessor(
    DataAccessorBase & data_accessor) {
  this->data_accessor = &data_accessor;
}

/* -------------------------------------------------------------------------- */
void SynchronizerRegistry::synchronize(SynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(data_accessor != nullptr, "No data accessor set.");

  std::pair<Tag2Sync::iterator, Tag2Sync::iterator> range =
      synchronizers.equal_range(tag);

  for (auto it = range.first; it != range.second; ++it) {
    it->second->synchronize(*data_accessor, tag);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SynchronizerRegistry::asynchronousSynchronize(SynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(data_accessor != nullptr, "No data accessor set.");

  std::pair<Tag2Sync::iterator, Tag2Sync::iterator> range =
      synchronizers.equal_range(tag);

  for (auto it = range.first; it != range.second; ++it) {
    (*it).second->asynchronousSynchronize(*data_accessor, tag);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SynchronizerRegistry::waitEndSynchronize(SynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(data_accessor != nullptr, "No data accessor set.");

  std::pair<Tag2Sync::iterator, Tag2Sync::iterator> range =
      synchronizers.equal_range(tag);

  for (auto it = range.first; it != range.second; ++it) {
    (*it).second->waitEndSynchronize(*data_accessor, tag);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SynchronizerRegistry::registerSynchronizer(Synchronizer & synchronizer,
                                                SynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  synchronizers.insert(
      std::pair<SynchronizationTag, Synchronizer *>(tag, &synchronizer));

  AKANTU_DEBUG_OUT();
}

} // namespace akantu
