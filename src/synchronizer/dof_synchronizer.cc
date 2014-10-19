/**
 * @file   dof_synchronizer.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
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
DOFSynchronizer::DOFSynchronizer(const Mesh & mesh, UInt nb_degree_of_freedom) :
  global_dof_equation_numbers(0, 1, "global_equation_number"),
  local_dof_equation_numbers(0, 1, "local_equation_number"),
  dof_global_ids(0, 1, "global_ids"),
  dof_types(0, 1, "types") {
  gather_scatter_scheme_initialized = false;

  communicator = &StaticCommunicator::getStaticCommunicator();

  prank = communicator->whoAmI();
  psize = communicator->getNbProc();

  proc_informations.resize(psize);

  UInt nb_nodes = mesh.getNbNodes();

  nb_dofs = nb_nodes * nb_degree_of_freedom;

  dof_global_ids.resize(nb_dofs);
  dof_types.resize(nb_dofs);

  nb_global_dofs = mesh.getNbGlobalNodes() * nb_degree_of_freedom;

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
    if(*dof_type >= 0) {
      proc_informations[*dof_type].master_dofs.push_back(n);
    }
    dof_type++;
  }

  /// exchanging information
  for (UInt p = 1; p < psize; ++p) {
    UInt sendto   = (psize + prank + p) % psize;
    UInt recvfrom = (psize + prank - p) % psize;

    UInt nb_master_dofs = proc_informations[sendto].master_dofs.getSize();
    UInt * master_dofs = proc_informations[sendto].master_dofs.storage();
    UInt * send_buffer = new UInt[nb_master_dofs];
    for (UInt d = 0; d < nb_master_dofs; ++d) {
      send_buffer[d] = dof_global_id[master_dofs[d]];
    }

    UInt nb_slave_dofs = 0;
    UInt * recv_buffer;
    std::vector<CommunicationRequest *> requests;
    requests.push_back(communicator->asyncSend(&nb_master_dofs, 1, sendto, 0));
    if(nb_master_dofs != 0) {
      AKANTU_DEBUG(dblInfo, "Sending "<< nb_master_dofs << " nodes to " << sendto + 1);
      requests.push_back(communicator->asyncSend(send_buffer, nb_master_dofs, sendto, 1));
    }
    communicator->receive(&nb_slave_dofs, 1, recvfrom, 0);
    if(nb_slave_dofs != 0) {
      AKANTU_DEBUG(dblInfo, "Receiving "<< nb_slave_dofs << " nodes from " << recvfrom + 1);
      proc_informations[recvfrom].slave_dofs.resize(nb_slave_dofs);
      recv_buffer = proc_informations[recvfrom].slave_dofs.storage();
      communicator->receive(recv_buffer, nb_slave_dofs, recvfrom, 1);
    }

    for (UInt d = 0; d < nb_slave_dofs; ++d) {
      recv_buffer[d] = global_dof_to_local[recv_buffer[d]];
    }

    communicator->waitAll(requests);
    communicator->freeCommunicationRequest(requests);
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
  communicator->allGather(nb_dof_per_proc, 1);
  AKANTU_DEBUG(dblDebug, "I have " << local_dofs->getSize() << " not ghost dofs ("
               << nb_local_dofs << " local and " << nb_needed_dofs << " slave)");

  UInt pos = 0;
  for (UInt p = 0; p < prank; ++p)  pos += nb_dof_per_proc[p];
  UInt nb_total_dofs = pos;
  for (UInt p = prank; p < psize; ++p)  nb_total_dofs += nb_dof_per_proc[p];

  Int * nb_values = new Int[psize];
  for (UInt p = 0; p < psize; ++p) nb_values[p] = nb_dof_per_proc[p] * 2;

  UInt * buffer = new UInt[2*nb_total_dofs];
  memcpy(buffer + 2 * pos, local_dofs->storage(), local_dofs->getSize() * 2 * sizeof(UInt));
  delete local_dofs;

  communicator->allGatherV(buffer, nb_values);

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
