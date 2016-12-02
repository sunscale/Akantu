/**
 * @file   synchronizer_tmpl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Wed Aug  3 13:49:36 2016
 *
 * @brief  Implementation of the helper classes for the synchronizer
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_array.hh"
#include "synchronizer.hh"
#include "synchronizer_impl.hh"
#include "data_accessor.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_SYNCHRONIZER_TMPL_HH__
#define __AKANTU_SYNCHRONIZER_TMPL_HH__

namespace akantu {

/// synchronize ghosts
template<class DataAccessorT>
void Synchronizer::synchronize(DataAccessorT & data_accessor,
                               const SynchronizationTag & tag) {
  if (SynchronizerImpl<Element> * synch_el =
          dynamic_cast<SynchronizerImpl<Element> *>(this)) {
    synch_el->synchronizeImpl(
        dynamic_cast< DataAccessor<Element> &>(data_accessor), tag);
  } else if (SynchronizerImpl<UInt> * synch_dof =
                 dynamic_cast<SynchronizerImpl<UInt> *>(this)) {
    synch_dof->synchronizeImpl(
        dynamic_cast< DataAccessor<UInt> &>(data_accessor), tag);
  } else {
    AKANTU_EXCEPTION("You synchronizer is not of a known type");
  }
}

/* -------------------------------------------------------------------------- */
template <class DataAccessorT>
void Synchronizer::asynchronousSynchronize(DataAccessorT & data_accessor,
                                           const SynchronizationTag & tag) {
  if (SynchronizerImpl<Element> * synch_el =
          dynamic_cast<SynchronizerImpl<Element> *>(this)) {
    synch_el->asynchronousSynchronizeImpl(
        dynamic_cast< DataAccessor<Element> &>(data_accessor), tag);
  } else if (SynchronizerImpl<UInt> * synch_dof =
                 dynamic_cast<SynchronizerImpl<UInt> *>(this)) {
    synch_dof->asynchronousSynchronizeImpl(
        dynamic_cast< DataAccessor<UInt> &>(data_accessor), tag);
  } else {
    AKANTU_EXCEPTION("You synchronizer is not of a known type");
  }
}

/* -------------------------------------------------------------------------- */
template <class DataAccessorT>
void Synchronizer::waitEndSynchronize(DataAccessorT & data_accessor,
                                      const SynchronizationTag & tag) {
  if (SynchronizerImpl<Element> * synch_el =
          dynamic_cast<SynchronizerImpl<Element> *>(this)) {
    synch_el->waitEndSynchronizeImpl(
        dynamic_cast< DataAccessor<Element> &>(data_accessor), tag);
  } else if (SynchronizerImpl<UInt> * synch_dof =
                 dynamic_cast<SynchronizerImpl<UInt> *>(this)) {
    synch_dof->waitEndSynchronizeImpl(
        dynamic_cast< DataAccessor<UInt> &>(data_accessor), tag);
  } else {
    AKANTU_EXCEPTION("You synchronizer is not of a known type");
  }
}

/// compute buffer size for a given tag and data accessor
template <class DataAccessorT>
void Synchronizer::computeBufferSize(DataAccessorT & data_accessor,
                                     const SynchronizationTag & tag) {
  if (SynchronizerImpl<Element> * synch_el =
          dynamic_cast<SynchronizerImpl<Element> *>(this)) {
    synch_el->computeBufferSizeImpl(
        dynamic_cast< DataAccessor<Element> &>(data_accessor), tag);
  } else if (SynchronizerImpl<UInt> * synch_dof =
                 dynamic_cast<SynchronizerImpl<UInt> *>(this)) {
    synch_dof->computeBufferSizeImpl(
        dynamic_cast< DataAccessor<UInt> &>(data_accessor), tag);
  } else {
    AKANTU_EXCEPTION("You synchronizer is not of a known type");
  }
}

} // akantu

#endif /* __AKANTU_SYNCHRONIZER_TMPL_HH__ */
