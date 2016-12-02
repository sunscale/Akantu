/**
 * @file   dof_synchronizer.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 17 2011
 * @date last modification: Wed Oct 21 2015
 *
 * @brief  DOF synchronizing object implementation
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
#include "dof_synchronizer.hh"
#include "dof_manager_default.hh"
#include "mesh.hh"
#include "node_synchronizer.hh"
/* -------------------------------------------------------------------------- */
#include <algorithm>
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
 * to be sent to prank and a second one with dofs that willl be received form
 * prank.
 * This information is needed for the asychronous communications. The
 * constructor sets up this information.
 */
DOFSynchronizer::DOFSynchronizer(DOFManagerDefault & dof_manager, const ID & id,
                                 MemoryID memory_id, StaticCommunicator & comm)
    : SynchronizerImpl<UInt>(id, memory_id, comm), root(0),
      dof_manager(dof_manager),
      slave_receive_dofs(0, 1, "dofs-to-receive-from-master"),
      dof_changed(true) {
  std::vector<ID> dof_ids = dof_manager.getDOFIDs();

  // Transfers nodes to global equation numbers in new schemes
  for (ID dof_id : dof_ids) {
    registerDOFs(dof_id);
  }

  this->initScatterGatherCommunicationScheme();
}

/* -------------------------------------------------------------------------- */
DOFSynchronizer::~DOFSynchronizer() {}

/* -------------------------------------------------------------------------- */
void DOFSynchronizer::registerDOFs(const ID & dof_id) {
  if(this->nb_proc == 1) return;

  typedef Communications<UInt>::const_scheme_iterator const_scheme_iterator;

  const Array<UInt> equation_numbers =
      dof_manager.getLocalEquationNumbers(dof_id);

  if (dof_manager.getSupportType(dof_id) == _dst_nodal) {
    const NodeSynchronizer & node_synchronizer =
        dof_manager.getMesh().getNodeSynchronizer();
    const Array<UInt> & associated_nodes =
        dof_manager.getDOFsAssociatedNodes(dof_id);

    const Communications<UInt> & node_communications =
        node_synchronizer.getCommunications();

    auto transcode_node_to_global_dof_scheme =
        [this, &associated_nodes,
         &equation_numbers](const_scheme_iterator it, const_scheme_iterator end,
                            const CommunicationSendRecv & sr) -> void {
      for (; it != end; ++it) {
        auto scheme = communications.createScheme(it->first, sr);

        auto node_scheme = it->second;
        for (auto node : node_scheme) {
          auto an_begin = associated_nodes.begin();
          auto an_it = an_begin;
          auto an_end = associated_nodes.begin();

          std::vector<UInt> global_dofs_per_node;
          while ((an_it = std::find(an_it, an_end, node)) != an_end) {
            UInt pos = an_it - an_end;
            UInt local_eq_num = equation_numbers(pos);
            UInt global_eq_num = local_eq_num;
            dof_manager.localToGlobalEquationNumber(global_eq_num);
            global_dofs_per_node.push_back(global_eq_num);
          }

          std::sort(global_dofs_per_node.begin(), global_dofs_per_node.end());
          std::transform(
              global_dofs_per_node.begin(), global_dofs_per_node.end(),
              global_dofs_per_node.begin(), [this](UInt g) -> UInt {
                UInt l = dof_manager.globalToLocalEquationNumber(g);
                return l;
              });
          for (auto leqnum : global_dofs_per_node) {
            scheme.push_back(leqnum);
          }
        }
      }
    };

    const_scheme_iterator nss_it = node_communications.begin_send_scheme();
    const_scheme_iterator nss_end = node_communications.end_send_scheme();

    transcode_node_to_global_dof_scheme(nss_it, nss_end, _send);

    const_scheme_iterator nrs_it = node_communications.begin_recv_scheme();
    const_scheme_iterator nrs_end = node_communications.end_recv_scheme();

    transcode_node_to_global_dof_scheme(nrs_it, nrs_end, _recv);
  }

  dof_changed = true;
}

/* -------------------------------------------------------------------------- */
void DOFSynchronizer::initScatterGatherCommunicationScheme() {
  AKANTU_DEBUG_IN();

  if (this->nb_proc == 1) {
    AKANTU_DEBUG_OUT();
    return;
  }

  UInt nb_dofs = dof_manager.getLocalSystemSize();

  this->slave_receive_dofs.clear();
  this->master_receive_dofs.clear();

  Array<UInt> dofs_to_send;
  for (UInt n = 0; n < nb_dofs; ++n) {
    if (dof_manager.isLocalOrMasterDOF(n)) {
      this->slave_receive_dofs.push_back(n);

      UInt global_dof = n;
      dof_manager.localToGlobalEquationNumber(global_dof);
      dofs_to_send.push_back(global_dof);
    }
  }

  if (this->rank == UInt(this->root)) {
    Array<UInt> nb_dof_per_proc(this->nb_proc);
    communicator.gather(dofs_to_send.getSize(), nb_dof_per_proc);

    std::vector<CommunicationRequest> requests;
    for (UInt p = 0; p < nb_proc; ++p) {
      if (p == UInt(this->root)) {
        continue;
      }
      Array<UInt> & receive_per_proc = master_receive_dofs[nb_proc];
      receive_per_proc.resize(nb_dof_per_proc(p));
      requests.push_back(communicator.asyncReceive(
          receive_per_proc, p,
          Tag::genTag(p, 0, Tag::_GATHER_INITIALIZATION, this->hash_id)));
    }

    communicator.waitAll(requests);
    communicator.freeCommunicationRequest(requests);
  } else {
    communicator.gather(dofs_to_send.getSize(), 0);
    AKANTU_DEBUG(dblDebug, "I have " << nb_dofs << " dofs ("
                                     << dofs_to_send.getSize()
                                     << " to send to master proc");

    communicator.send(
        dofs_to_send, 0,
        Tag::genTag(this->rank, 0, Tag::_GATHER_INITIALIZATION, this->hash_id));
  }

  dof_changed = false;

  AKANTU_DEBUG_OUT();
}

__END_AKANTU__
