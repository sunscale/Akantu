/**
 * @file   dof_synchronizer_inline_impl.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 17 2011
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  DOFSynchronizer inline implementation
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
#include "communication_buffer.hh"
#include "data_accessor.hh"
#include "dof_manager_default.hh"
#include "dof_synchronizer.hh"
/* -------------------------------------------------------------------------- */
#include <map>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_DOF_SYNCHRONIZER_INLINE_IMPL_CC__
#define __AKANTU_DOF_SYNCHRONIZER_INLINE_IMPL_CC__

namespace akantu {

/* -------------------------------------------------------------------------- */
template <typename T>
void DOFSynchronizer::gather(const Array<T> & to_gather, Array<T> & gathered) {
  if (this->hasChanged())
    initScatterGatherCommunicationScheme();

  AKANTU_DEBUG_ASSERT(this->rank == UInt(this->root),
                      "This function cannot be called on a slave processor");
  AKANTU_DEBUG_ASSERT(to_gather.size() ==
                          this->dof_manager.getLocalSystemSize(),
                      "The array to gather does not have the correct size");
  AKANTU_DEBUG_ASSERT(gathered.size() == this->dof_manager.getSystemSize(),
                      "The gathered array does not have the correct size");

  if (this->nb_proc == 1) {
    gathered.copy(to_gather, true);

    AKANTU_DEBUG_OUT();
    return;
  }

  std::map<UInt, CommunicationBuffer> buffers;
  std::vector<CommunicationRequest> requests;
  for (UInt p = 0; p < this->nb_proc; ++p) {
    if (p == UInt(this->root))
      continue;

    auto receive_it = this->master_receive_dofs.find(p);
    AKANTU_DEBUG_ASSERT(receive_it != this->master_receive_dofs.end(),
                        "Could not find the receive list for dofs of proc "
                            << p);
    const Array<UInt> & receive_dofs = receive_it->second;
    if (receive_dofs.size() == 0)
      continue;

    CommunicationBuffer & buffer = buffers[p];

    buffer.resize(receive_dofs.size() * to_gather.getNbComponent() * sizeof(T));

    AKANTU_DEBUG_INFO(
        "Preparing to receive data for "
        << receive_dofs.size() << " dofs from processor " << p << " "
        << Tag::genTag(p, this->root, Tag::_GATHER, this->hash_id));

    requests.push_back(communicator.asyncReceive(
        buffer, p, Tag::genTag(p, this->root, Tag::_GATHER, this->hash_id)));
  }

  auto data_gathered_it = gathered.begin(to_gather.getNbComponent());

  { // copy master data
    auto data_to_gather_it = to_gather.begin(to_gather.getNbComponent());
    for (auto local_dof : root_dofs) {
      UInt global_dof = dof_manager.localToGlobalEquationNumber(local_dof);

      Vector<T> dof_data_gathered = data_gathered_it[global_dof];
      Vector<T> dof_data_to_gather = data_to_gather_it[local_dof];
      dof_data_gathered = dof_data_to_gather;
    }
  }

  auto rr = UInt(-1);
  while ((rr = communicator.waitAny(requests)) != UInt(-1)) {
    CommunicationRequest & request = requests[rr];
    UInt sender = request.getSource();

    AKANTU_DEBUG_ASSERT(this->master_receive_dofs.find(sender) !=
                                this->master_receive_dofs.end() &&
                            buffers.find(sender) != buffers.end(),
                        "Missing infos concerning proc " << sender);

    const Array<UInt> & receive_dofs =
        this->master_receive_dofs.find(sender)->second;
    CommunicationBuffer & buffer = buffers[sender];

    for (auto global_dof : receive_dofs) {
      Vector<T> dof_data = data_gathered_it[global_dof];
      buffer >> dof_data;
    }

    requests.erase(requests.begin() + rr);
  }
}

/* -------------------------------------------------------------------------- */
template <typename T> void DOFSynchronizer::gather(const Array<T> & to_gather) {
  AKANTU_DEBUG_IN();

  if (this->hasChanged())
    initScatterGatherCommunicationScheme();

  AKANTU_DEBUG_ASSERT(this->rank != UInt(this->root),
                      "This function cannot be called on the root processor");
  AKANTU_DEBUG_ASSERT(to_gather.size() ==
                          this->dof_manager.getLocalSystemSize(),
                      "The array to gather does not have the correct size");

  if (this->root_dofs.size() == 0) {
    AKANTU_DEBUG_OUT();
    return;
  }
  CommunicationBuffer buffer(this->root_dofs.size() *
                             to_gather.getNbComponent() * sizeof(T));

  auto data_it = to_gather.begin(to_gather.getNbComponent());
  for (auto dof : this->root_dofs) {
    Vector<T> data = data_it[dof];
    buffer << data;
  }

  AKANTU_DEBUG_INFO("Gathering data for "
                    << to_gather.size() << " dofs on processor " << this->root
                    << " "
                    << Tag::genTag(this->rank, 0, Tag::_GATHER, this->hash_id));

  communicator.send(buffer, this->root,
                    Tag::genTag(this->rank, 0, Tag::_GATHER, this->hash_id));

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <typename T>
void DOFSynchronizer::scatter(Array<T> & scattered,
                              const Array<T> & to_scatter) {
  AKANTU_DEBUG_IN();

  if (this->hasChanged())
    initScatterGatherCommunicationScheme();
  AKANTU_DEBUG_ASSERT(this->rank == UInt(this->root),
                      "This function cannot be called on a slave processor");
  AKANTU_DEBUG_ASSERT(scattered.size() ==
                          this->dof_manager.getLocalSystemSize(),
                      "The scattered array does not have the correct size");
  AKANTU_DEBUG_ASSERT(to_scatter.size() == this->dof_manager.getSystemSize(),
                      "The array to scatter does not have the correct size");

  if (this->nb_proc == 1) {
    scattered.copy(to_scatter, true);
    AKANTU_DEBUG_OUT();
    return;
  }

  std::map<UInt, CommunicationBuffer> buffers;
  std::vector<CommunicationRequest> requests;

  for (UInt p = 0; p < nb_proc; ++p) {
    auto data_to_scatter_it = to_scatter.begin(to_scatter.getNbComponent());

    if (p == this->rank) {
      auto data_scattered_it = scattered.begin(to_scatter.getNbComponent());

      // copy the data for the local processor
      for (auto local_dof : root_dofs) {
        auto global_dof = dof_manager.localToGlobalEquationNumber(local_dof);

        Vector<T> dof_data_to_scatter = data_to_scatter_it[global_dof];
        Vector<T> dof_data_scattered = data_scattered_it[local_dof];
        dof_data_scattered = dof_data_to_scatter;
      }

      continue;
    }

    const Array<UInt> & receive_dofs =
        this->master_receive_dofs.find(p)->second;

    // prepare the send buffer
    CommunicationBuffer & buffer = buffers[p];
    buffer.resize(receive_dofs.size() * scattered.getNbComponent() * sizeof(T));

    // pack the data
    for (auto global_dof : receive_dofs) {
      Vector<T> dof_data_to_scatter = data_to_scatter_it[global_dof];
      buffer << dof_data_to_scatter;
    }

    // send the data
    requests.push_back(communicator.asyncSend(
        buffer, p, Tag::genTag(p, 0, Tag::_SCATTER, this->hash_id)));
  }

  // wait a clean communications
  communicator.waitAll(requests);
  communicator.freeCommunicationRequest(requests);

  // synchronize slave and ghost nodes
  synchronize(scattered);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <typename T> void DOFSynchronizer::scatter(Array<T> & scattered) {
  AKANTU_DEBUG_IN();

  if (this->hasChanged())
    this->initScatterGatherCommunicationScheme();
  AKANTU_DEBUG_ASSERT(this->rank != UInt(this->root),
                      "This function cannot be called on the root processor");
  AKANTU_DEBUG_ASSERT(scattered.size() ==
                          this->dof_manager.getLocalSystemSize(),
                      "The scattered array does not have the correct size");

  // prepare the data
  auto data_scattered_it = scattered.begin(scattered.getNbComponent());
  CommunicationBuffer buffer(this->root_dofs.size() *
                             scattered.getNbComponent() * sizeof(T));

  // receive the data
  communicator.receive(
      buffer, this->root,
      Tag::genTag(this->rank, 0, Tag::_SCATTER, this->hash_id));

  // unpack the data
  for (auto local_dof : root_dofs) {
    Vector<T> dof_data_scattered = data_scattered_it[local_dof];
    buffer >> dof_data_scattered;
  }

  // synchronize the ghosts
  synchronize(scattered);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <template <class> class Op, typename T>
void DOFSynchronizer::reduceSynchronize(Array<T> & array) const {
  ReduceDataAccessor<UInt, Op, T> data_accessor(array,
                                                SynchronizationTag::_whatever);
  this->slaveReductionOnceImpl(data_accessor, SynchronizationTag::_whatever);
  this->synchronize(array);
}

/* -------------------------------------------------------------------------- */
template <typename T>
void DOFSynchronizer::synchronize(Array<T> & array) const {
  SimpleUIntDataAccessor<T> data_accessor(array, SynchronizationTag::_whatever);
  this->synchronizeOnce(data_accessor, SynchronizationTag::_whatever);
}

} // namespace akantu

#endif /* __AKANTU_DOF_SYNCHRONIZER_INLINE_IMPL_CC__ */
