/**
 * @file   dof_synchronizer_inline_impl.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 17 2011
 * @date last modification: Tue Dec 09 2014
 *
 * @brief  DOFSynchronizer inline implementation
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include "communication_buffer.hh"
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
void DOFSynchronizer::gather(const Array<T> & to_gather,
                             Array<T> & gathered) {
  if (dof_changed)
    initScatterGatherCommunicationScheme();

  AKANTU_DEBUG_ASSERT(this->rank == UInt(this->root),
                      "This function cannot be called on a slave processor");
  AKANTU_DEBUG_ASSERT(to_gather.getSize() ==
                          this->dof_manager.getLocalSystemSize(),
                      "The array to gather does not have the correct size");
  AKANTU_DEBUG_ASSERT(gathered.getSize() == this->dof_manager.getSystemSize(),
                      "The gathered array does not have the correct size");

  if (this->nb_proc == 1) {
    gathered.copy(to_gather, true);

    AKANTU_DEBUG_OUT();
    return;
  }

  std::map<UInt, CommunicationBuffer> buffers;
  std::vector<CommunicationRequest> requests;
  for (UInt p = 1; p < nb_proc; ++p) {
    CommunicationBuffer & buffer = buffers[nb_proc];
    const Array<UInt> & receive_dofs =
        this->master_receive_dofs.find(p)->second;
    buffer.resize(receive_dofs.getSize() * to_gather.getNbComponent() *
                  sizeof(T));
    requests.push_back(communicator.asyncReceive(
        buffer, p, Tag::genTag(p, 0, Tag::_GATHER, this->hash_id)));
  }

  for (UInt p = 0; p < nb_proc; ++p) {
    auto data_it = to_gather.begin(to_gather.getNbComponent());

    if (p == this->rank) {
      auto gdata_it = gathered.begin(to_gather.getNbComponent());
      auto it = slave_receive_dofs.begin();
      auto end = slave_receive_dofs.end();
      for (; it != end; ++it) {
        Vector<T> svect = data_it[*it];

        UInt global_dof = *it;
        dof_manager.localToGlobalEquationNumber(global_dof);
        Vector<T> dvect = gdata_it[global_dof];
        dvect = svect;
      }
      continue;
    }
    UInt rr = communicator.waitAny(requests);

    CommunicationRequest & request = requests[rr];
    UInt sender = request.getSource();

    const Array<UInt> & receive_dofs =
        this->master_receive_dofs.find(sender)->second;
    CommunicationBuffer & buffer = buffers[sender];

    auto it = receive_dofs.begin();
    auto end = receive_dofs.end();
    for (; it != end; ++it) {
      Vector<T> vect = data_it[*it];
      buffer >> vect;
    }
  }

  communicator.freeCommunicationRequest(requests);
}

/* -------------------------------------------------------------------------- */
template <typename T>
void DOFSynchronizer::gather(const Array<T> & to_gather) {
  AKANTU_DEBUG_IN();

  if (dof_changed)
    initScatterGatherCommunicationScheme();

  AKANTU_DEBUG_ASSERT(this->rank != UInt(this->root),
                      "This function cannot be called on the root processor");
  AKANTU_DEBUG_ASSERT(to_gather.getSize() ==
                          this->dof_manager.getLocalSystemSize(),
                      "The array to gather does not have the correct size");

  CommunicationBuffer buffer(this->slave_receive_dofs.getSize() *
                             to_gather.getNbComponent() * sizeof(T));

  auto data_it = to_gather.begin(to_gather.getNbComponent());
  auto it = this->slave_receive_dofs.begin();
  auto end = this->slave_receive_dofs.end();

  for (; it != end; ++it) {
    buffer << data_it[*it];
  }

  communicator.send(buffer, this->root,
                    Tag::genTag(this->rank, 0, Tag::_GATHER, this->hash_id));

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <typename T>
void DOFSynchronizer::scatter(Array<T> & scattered,
                              const Array<T> & to_scatter) {
  AKANTU_DEBUG_IN();

  if (dof_changed)
    initScatterGatherCommunicationScheme();
  AKANTU_DEBUG_ASSERT(this->rank == UInt(this->root),
                      "This function cannot be called on a slave processor");
  AKANTU_DEBUG_ASSERT(scattered.getSize() ==
                          this->dof_manager.getLocalSystemSize(),
                      "The scattered array does not have the correct size");
  AKANTU_DEBUG_ASSERT(to_scatter.getSize() == this->dof_manager.getSystemSize(),
                      "The array to scatter does not have the correct size");

  if (this->nb_proc == 1) {
    scattered.copy(*to_scatter, true);
    AKANTU_DEBUG_OUT();
    return;
  }

  std::map<UInt, CommunicationBuffer> buffers;
  std::vector<CommunicationRequest> requests;

  for (UInt p = 0; p < nb_proc; ++p) {
    auto data_it = to_scatter.begin(to_scatter.getNbComponent());

    if (p == this->rank) {
      auto gdata_it = scattered.begin(to_scatter.getNbComponent());
      auto it = slave_receive_dofs.begin();
      auto end = slave_receive_dofs.end();
      for (; it != end; ++it) {
        UInt global_dof = *it;
        dof_manager.localToGlobalEquationNumber(global_dof);

        Vector<T> svect = data_it[global_dof];
        Vector<T> dvect = gdata_it[*it];
        dvect = svect;
      }
    }

    UInt rr = communicator.waitAny(requests);

    CommunicationRequest & request = requests[rr];
    UInt sender = request.getSource();

    const Array<UInt> & receive_dofs =
        this->master_receive_dofs.find(sender)->second;
    CommunicationBuffer & buffer = buffers[sender];

    buffer.resize(receive_dofs.getSize() * scattered.getNbComponent() *
                  sizeof(T));

    auto it = receive_dofs.begin();
    auto end = receive_dofs.end();
    for (; it != end; ++it) {
      Vector<T> vect = data_it[*it];
      buffer << vect;

      requests.push_back(communicator.asyncSend(
          buffer, p, Tag::genTag(p, 0, Tag::_SCATTER, this->hash_id)));
    }
  }

  communicator.freeCommunicationRequest(requests);

  // \todo copy local data

  synchronize(scattered);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <typename T>
void DOFSynchronizer::scatter(Array<T> & scattered) {
  if (dof_changed)
    this->initScatterGatherCommunicationScheme();
  AKANTU_DEBUG_ASSERT(this->rank != UInt(this->root),
                      "This function cannot be called on the root processor");
  AKANTU_DEBUG_ASSERT(scattered.getSize() ==
                          this->dof_manager.getLocalSystemSize(),
                      "The scattered array does not have the correct size");

  CommunicationBuffer buffer(this->slave_receive_dofs.getSize() *
                             scattered.getNbComponent() * sizeof(T));

  communicator.receive(
      buffer, this->root,
      Tag::genTag(this->rank, 0, Tag::_SCATTER, this->hash_id));

  auto data_it = scattered.begin(scattered.getNbComponent());
  auto it = this->slave_receive_dofs.begin();
  auto end = this->slave_receive_dofs.end();

  for (; it != end; ++it) {
    Vector<T> vect = data_it[*it];
    buffer >> vect;
  }

  synchronize(scattered);
  // synchronize the ghosts
}

/* -------------------------------------------------------------------------- */
template <typename T>
void DOFSynchronizer::synchronize(Array<T> & dof_vector) const {
  AKANTU_DEBUG_IN();

  if (this->nb_proc == 1) {
    AKANTU_DEBUG_OUT();
    return;
  }

  std::vector<CommunicationRequest> send_requests, recv_requests;
  std::vector<CommunicationBuffer> send_buffers, recv_buffers;

  auto rs_it = this->communications.begin_recv_scheme();
  auto rs_end = this->communications.end_recv_scheme();

  for (; rs_it != rs_end; ++rs_it) {
    auto scheme = rs_it->second;
    auto proc = rs_it->first;

    recv_buffers.emplace_back(
        {scheme.getSize() * dof_vector.getNbComponent() * sizeof(T)});
    recv_requests.push_back(communicator.asyncReceive(
        recv_buffers.back(), proc,
        Tag::genTag(this->rank, 0, Tag::_SYNCHRONIZE, this->hash_id)));
  }

  auto ss_it = this->communications.begin_send_scheme();
  auto ss_end = this->communications.end_send_scheme();

  for (; ss_it != ss_end; ++ss_it) {
    auto scheme = ss_it->second;
    auto proc = rs_it->first;

    send_buffers.emplace_back(
        {scheme.getSize() * dof_vector.getNbComponent() * sizeof(T)});
    CommunicationBuffer & buffer = send_buffers.back();

    auto data_it = dof_vector.begin(dof_vector.getNbComponent());
    for (auto dof : scheme) {
      buffer << data_it[dof];
    }

    send_requests.push_back(communicator.asyncSend(
        buffer, proc, Tag::genTag(proc, 0, Tag::_SYNCHRONIZE, this->hash_id)));
  }

  UInt ready;
  while ((ready = communicator.waitAny(recv_requests)) != UInt(-1)) {
    CommunicationRequest & req = recv_requests[ready];
    CommunicationBuffer & buffer = recv_buffers[ready];

    UInt proc = req.getSource();
    auto scheme = this->communications.getRecvScheme(proc);
    auto data_it = dof_vector.begin(dof_vector.getNbComponent());
    for (auto dof : scheme) {
      Vector<Real> vect = data_it[dof];
      buffer >> vect;
    }
  }

  communicator.waitAll(send_requests);
  communicator.freeCommunicationRequest(send_requests);
  communicator.freeCommunicationRequest(recv_requests);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <template <class> class Op, typename T>
void DOFSynchronizer::reduceSynchronize(Array<T> & dof_vector) const {
  AKANTU_DEBUG_IN();

  if (nb_proc == 1) {
    AKANTU_DEBUG_OUT();
    return;
  }

  std::vector<CommunicationRequest> send_requests, recv_requests;
  std::vector<CommunicationBuffer> send_buffers, recv_buffers;

  /// Same as synchronize but using send schemes as receive ones and vise versa
  auto rs_it = this->communications.begin_send_scheme();
  auto rs_end = this->communications.end_send_scheme();

  for (; rs_it != rs_end; ++rs_it) {
    auto scheme = rs_it->second;
    auto proc = rs_it->first;

    recv_buffers.emplace_back(
        {scheme.getSize() * dof_vector.getNbComponent() * sizeof(T)});
    recv_requests.push_back(communicator.asyncReceive(
        recv_buffers.back(), proc,
        Tag::genTag(this->rank, 0, Tag::_SYNCHRONIZE, this->hash_id)));
  }

  auto ss_it = this->communications.begin_recv_scheme();
  auto ss_end = this->communications.end_recv_scheme();

  for (; ss_it != ss_end; ++ss_it) {
    auto scheme = ss_it->second;
    auto proc = ss_it->first;

    send_buffers.emplace_back(
        {scheme.getSize() * dof_vector.getNbComponent() * sizeof(T)});
    CommunicationBuffer & buffer = send_buffers.back();

    auto data_it = dof_vector.begin(dof_vector.getNbComponent());
    for (auto dof : scheme) {
      buffer << data_it[dof];
    }

    send_requests.push_back(communicator.asyncSend(
        buffer, proc, Tag::genTag(proc, 0, Tag::_SYNCHRONIZE, this->hash_id)));
  }

  Op<Vector<T> > oper;

  UInt ready;
  while ((ready = communicator.waitAny(recv_requests)) != UInt(-1)) {
    CommunicationRequest & req = recv_requests[ready];
    CommunicationBuffer & buffer = recv_buffers[ready];

    UInt proc = req.getSource();
    auto scheme = this->communications.getRecvScheme(proc);
    auto data_it = dof_vector.begin(dof_vector.getNbComponent());

    Vector<Real> recv(dof_vector.getNbComponent());
    for (auto dof : scheme) {
      Vector<Real> vect = data_it[dof];
      buffer >> recv;

      oper(vect, recv);
    }
  }

  communicator.waitAll(send_requests);
  communicator.freeCommunicationRequest(send_requests);
  communicator.freeCommunicationRequest(recv_requests);

  synchronize(dof_vector);

  AKANTU_DEBUG_OUT();
}

} // akantu

#endif /* __AKANTU_DOF_SYNCHRONIZER_INLINE_IMPL_CC__ */
