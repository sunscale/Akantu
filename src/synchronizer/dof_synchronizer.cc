/**
 * @file   dof_synchronizer.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Aurelia Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 * @date creation: Fri Jun 17 2011
 * @date last modification: Thu Mar 27 2014
 *
 * @brief  DOF synchronizing object implementation
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
#include "dof_synchronizer.hh"
#include "mesh.hh"

/* -------------------------------------------------------------------------- */


__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
/**
 * A DOFSynchronizer needs a mesh and the number of degrees of freedom
 * per node to be created. In the constructor computes the local and global dof
 * number for each dof. The member
 * proc_informations (std vector) is resized with the number of mpi
 * processes. Each entry in the vector is a PerProcInformations object
 * that contains the interactions of the current mpi process (prank) with the
 * mpi process corresponding to the position of that entry. Every
 * ProcInformations object contains one array with the dofs that have
 * to be sent to prank and a second one with dofs that willl be received form prank.
 * This information is needed for the asychronous communications. The
 * constructor sets up this information.   
 * @param mesh mesh discretizing the domain we want to analyze 
 * @param nb_degree_of_freedom number of degrees of freedom per node
 */
DOFSynchronizer::DOFSynchronizer(const Mesh & mesh, UInt nb_degree_of_freedom) :
  global_dof_equation_numbers(0, 1, "global_equation_number"),
  local_dof_equation_numbers(0, 1, "local_equation_number"),
  dof_global_ids(0, 1, "global_ids"),
  dof_types(0, 1, "types") {

  gather_scatter_scheme_initialized = false;

  prank = static_communicator->whoAmI();
  psize = static_communicator->getNbProc();

  proc_informations.resize(psize);

  UInt nb_nodes = mesh.getNbNodes();

  nb_dofs = nb_nodes * nb_degree_of_freedom;

  dof_global_ids.resize(nb_dofs);
  dof_types.resize(nb_dofs);

  nb_global_dofs = mesh.getNbGlobalNodes() * nb_degree_of_freedom;

  /// compute the global id for each dof and store the dof type (pure ghost, slave, master or local)
  UInt * dof_global_id       = dof_global_ids.storage();
  Int  * dof_type            = dof_types.storage();
  for(UInt n = 0, ld = 0; n < nb_nodes; ++n) {
    UInt node_global_id = mesh.getNodeGlobalId(n);
    UInt node_type      = mesh.getNodeType(n);
    for (UInt d = 0; d < nb_degree_of_freedom; ++d, ++ld) {
      *dof_global_id = node_global_id * nb_degree_of_freedom + d;
      global_dof_to_local[*dof_global_id] = ld;

      *(dof_type++)  = node_type;

      dof_global_id++;
    }
  }

  if(psize == 1) return;

  /// creating communication scheme
  dof_global_id  = dof_global_ids.storage();
  dof_type       = dof_types.storage();

  for (UInt n = 0; n < nb_dofs; ++n) {
    /// check if dof is slave. In that case the value stored in
    /// dof_type corresponds to the process that has the corresponding
    /// master dof
    if(*dof_type >= 0) {
      /// has to receive n from proc[*dof_type]
      proc_informations[*dof_type].master_dofs.push_back(n);
    }
    dof_type++;
  }

  /// at this point the master nodes in PerProcInfo are known
  /// exchanging information
  for (UInt p = 1; p < psize; ++p) {
    UInt sendto   = (psize + prank + p) % psize;
    UInt recvfrom = (psize + prank - p) % psize;

    /// send the master nodes
    UInt nb_master_dofs = proc_informations[sendto].master_dofs.getSize();
    UInt * master_dofs = proc_informations[sendto].master_dofs.storage();
    UInt * send_buffer = new UInt[nb_master_dofs];
    for (UInt d = 0; d < nb_master_dofs; ++d) {
      send_buffer[d] = dof_global_id[master_dofs[d]];
    }

    UInt nb_slave_dofs = 0;
    UInt * recv_buffer;
    std::vector<CommunicationRequest *> requests;
    requests.push_back(static_communicator->asyncSend(&nb_master_dofs, 1, sendto, 0));
    if(nb_master_dofs != 0) {
      AKANTU_DEBUG(dblInfo, "Sending "<< nb_master_dofs << " dofs to " << sendto + 1);
      requests.push_back(static_communicator->asyncSend(send_buffer, nb_master_dofs, sendto, 1));
    }

    /// Receive the info and store them as slave nodes 
    static_communicator->receive(&nb_slave_dofs, 1, recvfrom, 0);
    if(nb_slave_dofs != 0) {
      AKANTU_DEBUG(dblInfo, "Receiving "<< nb_slave_dofs << " dofs from " << recvfrom + 1);
      proc_informations[recvfrom].slave_dofs.resize(nb_slave_dofs);
      recv_buffer = proc_informations[recvfrom].slave_dofs.storage();
      static_communicator->receive(recv_buffer, nb_slave_dofs, recvfrom, 1);
    }

    for (UInt d = 0; d < nb_slave_dofs; ++d) {
      recv_buffer[d] = global_dof_to_local[recv_buffer[d]];
    }

    static_communicator->waitAll(requests);
    static_communicator->freeCommunicationRequest(requests);
    requests.clear();
    delete [] send_buffer;
  }
}

/* -------------------------------------------------------------------------- */
DOFSynchronizer::~DOFSynchronizer() {

}


/* -------------------------------------------------------------------------- */
void DOFSynchronizer::initLocalDOFEquationNumbers() {
  AKANTU_DEBUG_IN();
  local_dof_equation_numbers.resize(nb_dofs);

  Int  * dof_equation_number = local_dof_equation_numbers.storage();
  for (UInt d = 0; d < nb_dofs; ++d) {
    *(dof_equation_number++) = d;
  }

  //local_dof_equation_numbers.resize(nb_dofs);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DOFSynchronizer::asynchronousSynchronize(DataAccessor & data_accessor,
					      SynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  if (communications.find(tag) == communications.end())
    computeBufferSize(data_accessor, tag);

  Communication & communication = communications[tag];

  AKANTU_DEBUG_ASSERT(communication.send_requests.size() == 0,
		      "There must be some pending sending communications. Tag is " << tag);

  std::map<SynchronizationTag, UInt>::iterator t_it = tag_counter.find(tag);
  UInt counter = 0;
  if(t_it == tag_counter.end()) {
    tag_counter[tag] = 0;
  } else {
    counter = ++(t_it->second);
  }

  for (UInt p = 0; p < psize; ++p) {
    UInt ssize = communication.size_to_send[p];
    if(p == prank || ssize == 0) continue;

    CommunicationBuffer & buffer = communication.send_buffer[p];
    buffer.resize(ssize);
#ifndef AKANTU_NDEBUG
    UInt nb_dofs   =  proc_informations[p].slave_dofs.getSize();
    AKANTU_DEBUG_INFO("Packing data for proc " << p
		      << " (" << ssize << "/" << nb_dofs
		      <<" data to send/dofs)");

    /// pack global equation numbers in debug mode
    Array<UInt>::const_iterator<UInt> bit  = proc_informations[p].slave_dofs.begin();
    Array<UInt>::const_iterator<UInt> bend = proc_informations[p].slave_dofs.end();
    for (; bit != bend; ++bit) {
      buffer << global_dof_equation_numbers[*bit];
    }
#endif


    /// dof synchronizer needs to send the data corresponding to the 
    data_accessor.packDOFData(buffer, proc_informations[p].slave_dofs, tag);

    AKANTU_DEBUG_ASSERT(buffer.getPackedSize() == ssize,
			"a problem has been introduced with "
			<< "false sent sizes declaration "
			<< buffer.getPackedSize() << " != " << ssize);
    AKANTU_DEBUG_INFO("Posting send to proc " << p
		      << " (tag: " << tag << " - " << ssize << " data to send)"
		      << " [" << Tag::genTag(prank, counter, tag) << "]");
    communication.send_requests.push_back(static_communicator->asyncSend(buffer.storage(),
									 ssize,
									 p,
									 Tag::genTag(prank, counter, tag)));
  }

  AKANTU_DEBUG_ASSERT(communication.recv_requests.size() == 0,
		      "There must be some pending receive communications");

  for (UInt p = 0; p < psize; ++p) {
    UInt rsize = communication.size_to_receive[p];
    if(p == prank || rsize == 0) continue;
    CommunicationBuffer & buffer = communication.recv_buffer[p];
    buffer.resize(rsize);

    AKANTU_DEBUG_INFO("Posting receive from proc " << p
		      << " (tag: " << tag << " - " << rsize << " data to receive) "
		      << " [" << Tag::genTag(p, counter, tag) << "]");
    communication.recv_requests.push_back(static_communicator->asyncReceive(buffer.storage(),
									    rsize,
									    p,
									    Tag::genTag(p, counter, tag)));
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DOFSynchronizer::waitEndSynchronize(DataAccessor & data_accessor,
					 SynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(communications.find(tag) != communications.end(), "No communication with the tag \""
		      << tag <<"\" started");

  Communication & communication = communications[tag];

  std::vector<CommunicationRequest *> req_not_finished;
  std::vector<CommunicationRequest *> * req_not_finished_tmp = &req_not_finished;
  std::vector<CommunicationRequest *> * recv_requests_tmp = &(communication.recv_requests);

  //  static_communicator->waitAll(recv_requests);
  while(!recv_requests_tmp->empty()) {
    for (std::vector<CommunicationRequest *>::iterator req_it = recv_requests_tmp->begin();
	 req_it != recv_requests_tmp->end() ; ++req_it) {
      CommunicationRequest * req = *req_it;

      if(static_communicator->testRequest(req)) {
	UInt proc = req->getSource();
	AKANTU_DEBUG_INFO("Unpacking data coming from proc " << proc);
	CommunicationBuffer & buffer = communication.recv_buffer[proc];

#ifndef AKANTU_NDEBUG
	Array<UInt>::const_iterator<UInt> bit  = proc_informations[proc].master_dofs.begin();
	Array<UInt>::const_iterator<UInt> bend = proc_informations[proc].master_dofs.end();

	for (; bit != bend; ++bit) {
	  Int global_dof_eq_nb_loc = global_dof_equation_numbers[*bit];
	  Int global_dof_eq_nb = 0;
	  buffer >> global_dof_eq_nb;
	  Real tolerance = Math::getTolerance();
	  if(std::abs(global_dof_eq_nb - global_dof_eq_nb_loc) <= tolerance) continue;
	  AKANTU_DEBUG_ERROR("Unpacking an unknown global dof equation number for dof: "
			     << *bit
			     << "(global dof equation number = " << global_dof_eq_nb
			     << " and buffer = " << global_dof_eq_nb << ") ["
			     << std::abs(global_dof_eq_nb - global_dof_eq_nb_loc)
			     << "] - tag: " << tag);
	}
	
#endif

	data_accessor.unpackDOFData(buffer, proc_informations[proc].master_dofs, tag);
	buffer.resize(0);

	AKANTU_DEBUG_ASSERT(buffer.getLeftToUnpack() == 0,
			    "all data have not been unpacked: "
			    << buffer.getLeftToUnpack() << " bytes left");
	static_communicator->freeCommunicationRequest(req);
      } else {
	req_not_finished_tmp->push_back(req);
      }
    }

    std::vector<CommunicationRequest *> * swap = req_not_finished_tmp;
    req_not_finished_tmp = recv_requests_tmp;
    recv_requests_tmp = swap;

    req_not_finished_tmp->clear();
  }


  AKANTU_DEBUG_INFO("Waiting that every send requests are received");
  static_communicator->waitAll(communication.send_requests);
  for (std::vector<CommunicationRequest *>::iterator req_it = communication.send_requests.begin();
       req_it != communication.send_requests.end() ; ++req_it) {
    CommunicationRequest & req = *(*req_it);

    if(static_communicator->testRequest(&req)) {
      UInt proc = req.getDestination();
      CommunicationBuffer & buffer = communication.send_buffer[proc];
      buffer.resize(0);
      static_communicator->freeCommunicationRequest(&req);
    }
  }
  communication.send_requests.clear();

  AKANTU_DEBUG_OUT();				


}

/* -------------------------------------------------------------------------- */
/**
 * This function computes the buffer size needed for in order to send data or receive data. 
 * @param data_accessor object of a class that needs to use the DOFSynchronizer. This class has to inherit from the  * DataAccessor interface.
 * @param tag synchronization tag: indicates what variable should be sychronized, e.g. mass, nodal residual, etc. 
 * for the SolidMechanicsModel
 */
void DOFSynchronizer::computeBufferSize(DataAccessor & data_accessor,
					SynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  communications[tag].resize(psize);

  for (UInt p = 0; p < psize; ++p) {
    /// initialize the size of data that we be sent and received
    UInt ssend    = 0;
    UInt sreceive = 0;
    if(p != prank) {
      /// check if processor prank has to send dof information to p
      if(proc_informations[p].slave_dofs.getSize() != 0) {

       
#ifndef AKANTU_NDEBUG
	/// in debug mode increase buffer size to send the positions
	/// of nodes in the direction that correspond to the dofs
	ssend += proc_informations[p].slave_dofs.getSize() * sizeof(Int);
#endif
	ssend += data_accessor.getNbDataForDOFs(proc_informations[p].slave_dofs, tag);
	AKANTU_DEBUG_INFO("I have " << ssend << "(" << ssend / 1024.
			  << "kB - "<< proc_informations[p].slave_dofs.getSize() <<" dof(s)) data to send to " << p << " for tag "
			  << tag);
      }

      /// check if processor prank has to receive dof information from p
      if(proc_informations[p].master_dofs.getSize() != 0) {
#ifndef AKANTU_NDEBUG
	/// in debug mode increase buffer size to receive the
	/// positions of nodes in the direction that correspond to the
	/// dofs
	sreceive += proc_informations[p].master_dofs.getSize() * sizeof(Int);
#endif
	sreceive += data_accessor.getNbDataForDOFs(proc_informations[p].master_dofs, tag);
	AKANTU_DEBUG_INFO("I have " << sreceive << "(" << sreceive / 1024.
			  << "kB - "<< proc_informations[p].master_dofs.getSize() <<" dof(s)) data to receive for tag "
			  << tag);
      }
    }

    communications[tag].size_to_send   [p] = ssend;
    communications[tag].size_to_receive[p] = sreceive;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DOFSynchronizer::initGlobalDOFEquationNumbers() {
  AKANTU_DEBUG_IN();
  global_dof_equation_numbers.resize(nb_dofs);

  Int  * dof_type            = dof_types.storage();
  UInt * dof_global_id       = dof_global_ids.storage();
  Int  * dof_equation_number = global_dof_equation_numbers.storage();
  for(UInt d = 0; d < nb_dofs; ++d) {
    /// if ghost dof the equation_number is greater than nb_global_dofs
    Int global_eq_num = *dof_global_id + (*dof_type > -3 ? 0 : nb_global_dofs);
    *(dof_equation_number)  = global_eq_num;
    global_dof_equation_number_to_local[global_eq_num] = d;

    dof_equation_number++;
    dof_global_id++;
    dof_type++;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DOFSynchronizer::initScatterGatherCommunicationScheme() {
  AKANTU_DEBUG_IN();

  if(psize == 1 || gather_scatter_scheme_initialized) {
    if(psize == 1) gather_scatter_scheme_initialized = true;
    AKANTU_DEBUG_OUT();
    return;
  }

  /// creating communication scheme
  UInt * dof_global_id = dof_global_ids.storage();
  Int * dof_type       = dof_types.storage();

  Array<UInt> * local_dofs = new Array<UInt>(0,2);
  UInt local_dof_val[2];
  nb_needed_dofs = 0;
  nb_local_dofs = 0;

  for (UInt n = 0; n < nb_dofs; ++n) {
    if(*dof_type != -3) {
      local_dof_val[0] = *dof_global_id;

      if(*dof_type == -1 || *dof_type == -2) {
        local_dof_val[1] = 0; // master for this node, shared or not
        nb_local_dofs++;
      } else if (*dof_type >= 0) {
        nb_needed_dofs++;
        local_dof_val[1] = 1; // slave node
      }
      local_dofs->push_back(local_dof_val);
    }
    dof_type++;
    dof_global_id++;
  }

  Int * nb_dof_per_proc = new Int[psize];
  nb_dof_per_proc[prank] = local_dofs->getSize();
  static_communicator->allGather(nb_dof_per_proc, 1);
  AKANTU_DEBUG(dblDebug, "I have " << local_dofs->getSize() << " not ghost dofs ("
               << nb_local_dofs << " local and " << nb_needed_dofs << " slave)");

  UInt pos = 0;
  for (UInt p = 0; p < prank; ++p)  pos += nb_dof_per_proc[p];
  UInt nb_total_dofs = pos;
  for (UInt p = prank; p < psize; ++p)  nb_total_dofs += nb_dof_per_proc[p];

  int * nb_values = new int[psize];
  for (unsigned int p = 0; p < psize; ++p) nb_values[p] = nb_dof_per_proc[p] * 2;

  UInt * buffer = new UInt[2*nb_total_dofs];
  memcpy(buffer + 2 * pos, local_dofs->storage(), local_dofs->getSize() * 2 * sizeof(UInt));
  delete local_dofs;

  static_communicator->allGatherV(buffer, nb_values);

  UInt * tmp_buffer = buffer;

  for (UInt p = 0; p < psize; ++p) {
    UInt proc_p_nb_dof = nb_dof_per_proc[p];
    if (p != prank){
      AKANTU_DEBUG(dblDebug, "I get " << proc_p_nb_dof << "(" << nb_values[p] << ") dofs from " << p + 1);
      proc_informations[p].dofs.resize(0);
      for (UInt dd = 0; dd < proc_p_nb_dof; ++dd) {
        UInt dof = tmp_buffer[2*dd];
        if(tmp_buffer[2*dd + 1] == 0)
          proc_informations[p].dofs.push_back(dof);
        else
          proc_informations[p].needed_dofs.push_back(dof);
      }

      AKANTU_DEBUG(dblDebug, "Proc " << p + 1 << " sends me "
                   << proc_informations[p].dofs.getSize()	<< " local dofs, and "
                   << proc_informations[p].needed_dofs.getSize() << " slave dofs");
    }
    tmp_buffer += 2*proc_p_nb_dof;
  }
  delete [] nb_dof_per_proc;
  delete [] buffer;
  delete [] nb_values;
  gather_scatter_scheme_initialized = true;

  AKANTU_DEBUG_OUT();
}

__END_AKANTU__
