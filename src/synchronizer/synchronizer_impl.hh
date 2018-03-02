/**
 * @file   synchronizer_impl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Implementation of the generic part of synchronizers
 *
 * @section LICENSE
 *
 * Copyright (©) 2016-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "communications.hh"
#include "synchronizer.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_SYNCHRONIZER_IMPL_HH__
#define __AKANTU_SYNCHRONIZER_IMPL_HH__

namespace akantu {

template <class Entity> class SynchronizerImpl : public Synchronizer {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  SynchronizerImpl(const Communicator & communicator,
                   const ID & id = "synchronizer", MemoryID memory_id = 0);

  SynchronizerImpl(const SynchronizerImpl & other, const ID & id);

  ~SynchronizerImpl() override = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// synchronous synchronization without state
  virtual void slaveReductionOnceImpl(DataAccessor<Entity> & data_accessor,
                                      const SynchronizationTag & tag) const;

  /// synchronous synchronization without state
  virtual void synchronizeOnceImpl(DataAccessor<Entity> & data_accessor,
                                   const SynchronizationTag & tag) const;

  /// asynchronous synchronization of ghosts
  virtual void
  asynchronousSynchronizeImpl(const DataAccessor<Entity> & data_accessor,
                              const SynchronizationTag & tag);

  /// wait end of asynchronous synchronization of ghosts
  virtual void waitEndSynchronizeImpl(DataAccessor<Entity> & data_accessor,
                                      const SynchronizationTag & tag);

  /// compute all buffer sizes
  virtual void
  computeAllBufferSizes(const DataAccessor<Entity> & data_accessor);

  /// compute buffer size for a given tag and data accessor
  virtual void computeBufferSizeImpl(const DataAccessor<Entity> & data_accessor,
                                     const SynchronizationTag & tag);

  /* ------------------------------------------------------------------------ */
  virtual void synchronizeImpl(DataAccessor<Entity> & data_accessor,
                               const SynchronizationTag & tag) {
    this->asynchronousSynchronizeImpl(data_accessor, tag);
    this->waitEndSynchronizeImpl(data_accessor, tag);
  }

  /* ------------------------------------------------------------------------ */
  /// extract the elements that have a true predicate from in_synchronizer and
  /// store them in the current synchronizer
  template <typename Pred>
  void split(SynchronizerImpl & in_synchronizer, Pred && pred);

  /// update schemes in a synchronizer
  template <typename Updater> void updateSchemes(Updater && scheme_updater);

  /// flip send and receive schemes
  void swapSendRecv();

  /* ------------------------------------------------------------------------ */
  virtual UInt sanityCheckDataSize(const Array<Entity> & elements,
                                   const SynchronizationTag & tag) const;
  virtual void
  packSanityCheckData(CommunicationDescriptor<Entity> & comm_desc) const;
  virtual void
  unpackSanityCheckData(CommunicationDescriptor<Entity> & comm_desc) const;

public:
  AKANTU_GET_MACRO(Communications, communications,
                   const Communications<Entity> &);

protected:
  AKANTU_GET_MACRO_NOT_CONST(Communications, communications,
                             Communications<Entity> &);

  virtual Int getRank(const Entity & entity) const = 0;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// information on the communications
  Communications<Entity> communications;
};

} // namespace akantu

#include "synchronizer_impl_tmpl.hh"

#endif /* __AKANTU_SYNCHRONIZER_IMPL_HH__ */
