/**
 * @file   dof_synchronizer_inline_impl.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 17 2011
 * @date last modification: Thu Jun 05 2014
 *
 * @brief  DOFSynchronizer inline implementation
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


/* -------------------------------------------------------------------------- */
template<typename T> void DOFSynchronizer::gather(const Array<T> & to_gather, UInt root,
						  Array<T> * gathered) const {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(gather_scatter_scheme_initialized == true, "You should call initScatterGatherCommunicationScheme before trying to gather a Array !!");

  if(psize == 1) {
    gathered->copy(to_gather, true);

    AKANTU_DEBUG_OUT();
    return;
  }

  UInt * dof_global_id = dof_global_ids.storage();
  Int  * dof_type      = dof_types.storage();
  T * to_gather_val    = to_gather.storage();

  T  * buffer;
  if(prank == root) {
    gathered->resize(nb_global_dofs / gathered->getNbComponent());
    buffer = gathered->storage();
  } else
    buffer = new T[nb_local_dofs];

  UInt dof = 0;
  for (UInt d = 0; d < nb_dofs; ++d) {
    if(*dof_type != -3) {
      if(*dof_type == -1 || *dof_type == -2) {
	if (prank == root) dof = *dof_global_id;
	buffer[dof++] = *to_gather_val;
      }
    }
    dof_type++;
    dof_global_id++;
    to_gather_val++;
  }

  if (prank == root) {
    for (UInt p = 0; p < psize; ++p) {
      if(p == root) continue;
      UInt nb_dofs = proc_informations[p].dofs.getSize();
      AKANTU_DEBUG_INFO("Gather - Receiving " << nb_dofs << " from " << p);
      if(nb_dofs) {
	buffer = new T[nb_dofs];
	communicator->receive(buffer, nb_dofs, p, 0);

	T * buffer_tmp = buffer;
	UInt * remote_dofs = proc_informations[p].dofs.storage();
	for (UInt d = 0; d < nb_dofs; ++d) {
	  gathered->storage()[*remote_dofs++] = *(buffer_tmp++);
	}

	delete [] buffer;
      }
    }
  } else {
    AKANTU_DEBUG_INFO("Gather - Sending " << nb_local_dofs << " to " << root);
    if(nb_local_dofs)
      communicator->send(buffer, nb_local_dofs, root, 0);
    delete [] buffer;
  }



  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<typename T> void DOFSynchronizer::scatter(Array<T> & scattered, UInt root,
				  const Array<T> * to_scatter) const {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(gather_scatter_scheme_initialized == true, "You should call initScatterGatherCommunicationScheme before trying to scatter a Array !!");

  if(psize == 1) {
    scattered.copy(*to_scatter, true);
    AKANTU_DEBUG_OUT();
    return;
  }

  scattered.resize(nb_dofs / scattered.getNbComponent());

  UInt * dof_global_id = dof_global_ids.storage();
  Int  * dof_type      = dof_types.storage();

  if (prank == root) {
    for (UInt p = 0; p < psize; ++p) {
      if(p == root) continue;

      UInt nb_needed_dof = proc_informations[p].needed_dofs.getSize();
      UInt nb_local_dof  = proc_informations[p].dofs.getSize();

      AKANTU_DEBUG_INFO("Scatter - Sending " << nb_local_dof + nb_needed_dof << " to " << p);
      if(nb_local_dof + nb_needed_dof) {
	T * send_buffer = new T[nb_local_dof + nb_needed_dof];

	T * buffer_tmp = send_buffer;

	UInt * remote_dofs = proc_informations[p].dofs.storage();
	for (UInt d = 0; d < nb_local_dof; ++d) {
	  *(buffer_tmp++) = to_scatter->storage()[*remote_dofs++];
	}

	remote_dofs = proc_informations[p].needed_dofs.storage();
	for (UInt d = 0; d < nb_needed_dof; ++d) {
	  *(buffer_tmp++) = to_scatter->storage()[*remote_dofs++];
	}

	communicator->send(send_buffer, nb_local_dof + nb_needed_dof, p, 0);
	delete [] send_buffer;
      }
    }

    T * scattered_val = scattered.storage();
    for (UInt d = 0; d < nb_dofs; ++d) {
      if(*dof_type != -3) {
	if(*dof_type >= -2)
	  *scattered_val = to_scatter->storage()[*dof_global_id];
      }
      scattered_val++;
      dof_type++;
      dof_global_id++;
    }
  } else {
    T  * buffer;
    AKANTU_DEBUG_INFO("Scatter - Receiving " << nb_dofs << " from " << root);
    if(nb_local_dofs + nb_needed_dofs) {
      buffer = new T[nb_local_dofs + nb_needed_dofs];
      communicator->receive(buffer, nb_local_dofs + nb_needed_dofs, root, 0);


      T * scattered_val = scattered.storage();
      UInt local_dofs = 0;
      UInt needed_dofs = nb_local_dofs;
      for (UInt d = 0; d < nb_dofs; ++d) {
	if(*dof_type != -3) {
	  if(*dof_type == -1 || *dof_type == -2)
	    *scattered_val = buffer[local_dofs++];
	  else if(*dof_type >= 0)
	    *scattered_val = buffer[needed_dofs++];
	}
	scattered_val++;
	dof_type++;
      }
      delete [] buffer;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<typename T> void DOFSynchronizer::synchronize(Array<T> & dof_vector) const {
  AKANTU_DEBUG_IN();

  if(psize == 1) {
    AKANTU_DEBUG_OUT();
    return;
  }


  for (UInt p = 1; p < psize; ++p) {
    UInt sendto   = (psize + prank + p) % psize;
    UInt recvfrom = (psize + prank - p) % psize;

    UInt nb_slave_dofs = proc_informations[sendto].slave_dofs.getSize();
    UInt * slave_dofs  = proc_informations[sendto].slave_dofs.storage();

    UInt nb_master_dofs = proc_informations[recvfrom].master_dofs.getSize();
    UInt * master_dofs  = proc_informations[recvfrom].master_dofs.storage();

    T * send_buffer = new T[nb_slave_dofs];
    T * recv_buffer = new T[nb_master_dofs];

    for (UInt d = 0; d < nb_slave_dofs; ++d) {
      AKANTU_DEBUG_ASSERT(dof_types(slave_dofs[d]) == -2,
			  "Sending node " << slave_dofs[d]
			  << "(gid" << dof_global_ids(d) << ") to proc "
			  << sendto << " but it is not a master node.");
      send_buffer[d] = dof_vector.storage()[slave_dofs[d]];
    }

    /// ring blocking communications
    CommunicationRequest * request = NULL;
    if(nb_slave_dofs  != 0) request = communicator->asyncSend(send_buffer, nb_slave_dofs,  sendto  , 0);
    if(nb_master_dofs != 0) communicator->receive(recv_buffer, nb_master_dofs, recvfrom, 0);

    for (UInt d = 0; d < nb_master_dofs; ++d) {
      AKANTU_DEBUG_ASSERT(dof_types(master_dofs[d]) >= 0,
			  "Received node " << master_dofs[d]
			  << "(gid" << dof_global_ids(d) << ")  from proc "
			  << recvfrom << " but it is not a slave node.");
      dof_vector.storage()[master_dofs[d]] = recv_buffer[d];
    }

    if(nb_slave_dofs != 0) {
      communicator->wait(request);
      communicator->freeCommunicationRequest(request);
    }
    delete [] send_buffer;
    delete [] recv_buffer;
  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
template<template <class> class Op, typename T> void DOFSynchronizer::reduceSynchronize(Array<T> & dof_vector) const {
  AKANTU_DEBUG_IN();

  if(psize == 1) {
    AKANTU_DEBUG_OUT();
    return;
  }

  for (UInt p = 1; p < psize; ++p) {
    UInt sendto   = (psize + prank + p) % psize;
    UInt recvfrom = (psize + prank - p) % psize;

    UInt nb_slave_dofs = proc_informations[recvfrom].slave_dofs.getSize();
    UInt * slave_dofs  = proc_informations[recvfrom].slave_dofs.storage();

    UInt nb_master_dofs = proc_informations[sendto].master_dofs.getSize();
    UInt * master_dofs  = proc_informations[sendto].master_dofs.storage();

    T * send_buffer = new T[nb_master_dofs];
    T * recv_buffer = new T[nb_slave_dofs];

    for (UInt d = 0; d < nb_master_dofs; ++d) {
      send_buffer[d] = dof_vector.storage()[master_dofs[d]];
    }

    CommunicationRequest * request = NULL;
    if(nb_master_dofs != 0) request = communicator->asyncSend(send_buffer, nb_master_dofs, sendto  , 0);
    if(nb_slave_dofs  != 0) communicator->receive(recv_buffer, nb_slave_dofs,  recvfrom, 0);

    Op<T> oper;

    for (UInt d = 0; d < nb_slave_dofs; ++d) {
      dof_vector.storage()[slave_dofs[d]] = oper(dof_vector.storage()[slave_dofs[d]], recv_buffer[d]);
    }

    if(nb_master_dofs != 0) {
      communicator->wait(request);
      communicator->freeCommunicationRequest(request);
    }
    delete [] send_buffer;
    delete [] recv_buffer;
  }

  synchronize(dof_vector);

  AKANTU_DEBUG_OUT();
}
