/**
 * @file   dof_synchronizer.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 17 2011
 * @date last modification: Tue Feb 06 2018
 *
 * @brief  DOF synchronizing object implementation
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
#include "dof_synchronizer.hh"
#include "aka_iterators.hh"
#include "dof_manager_default.hh"
#include "mesh.hh"
#include "node_synchronizer.hh"
/* -------------------------------------------------------------------------- */
#include <algorithm>
/* -------------------------------------------------------------------------- */

namespace akantu {

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
                                 MemoryID memory_id)
    : SynchronizerImpl<UInt>(dof_manager.getCommunicator(), id, memory_id),
      dof_manager(dof_manager) {
  std::vector<ID> dof_ids = dof_manager.getDOFIDs();

  // Transfers nodes to global equation numbers in new schemes
  for (const ID & dof_id : dof_ids) {
    registerDOFs(dof_id);
  }
}

/* -------------------------------------------------------------------------- */
DOFSynchronizer::~DOFSynchronizer() = default;

/* -------------------------------------------------------------------------- */
void DOFSynchronizer::registerDOFs(const ID & dof_id) {
  if (this->nb_proc == 1)
    return;

  if (dof_manager.getSupportType(dof_id) != _dst_nodal)
    return;

  const auto & equation_numbers = dof_manager.getLocalEquationsNumbers(dof_id);

  const auto & associated_nodes = dof_manager.getDOFsAssociatedNodes(dof_id);
  const auto & node_synchronizer = dof_manager.getMesh().getNodeSynchronizer();
  const auto & node_communications = node_synchronizer.getCommunications();

  auto transcode_node_to_global_dof_scheme =
      [this, &associated_nodes, &equation_numbers](
          auto && it, auto && end, const CommunicationSendRecv & sr) -> void {
    for (; it != end; ++it) {
      auto & scheme = communications.createScheme(it->first, sr);

      const auto & node_scheme = it->second;
      for (auto & node : node_scheme) {
        auto an_begin = associated_nodes.begin();
        auto an_it = an_begin;
        auto an_end = associated_nodes.end();

        std::vector<UInt> global_dofs_per_node;
        while ((an_it = std::find(an_it, an_end, node)) != an_end) {
          UInt pos = an_it - an_begin;
          UInt local_eq_num = equation_numbers(pos);
          UInt global_eq_num =
              dof_manager.localToGlobalEquationNumber(local_eq_num);
          global_dofs_per_node.push_back(global_eq_num);
          ++an_it;
        }

        std::sort(global_dofs_per_node.begin(), global_dofs_per_node.end());
        std::transform(global_dofs_per_node.begin(), global_dofs_per_node.end(),
                       global_dofs_per_node.begin(), [this](UInt g) -> UInt {
                         UInt l = dof_manager.globalToLocalEquationNumber(g);
                         return l;
                       });
        for (auto & leqnum : global_dofs_per_node) {
          scheme.push_back(leqnum);
        }
      }
    }
  };

  for (auto sr_it = send_recv_t::begin(); sr_it != send_recv_t::end();
       ++sr_it) {
    auto ncs_it = node_communications.begin_scheme(*sr_it);
    auto ncs_end = node_communications.end_scheme(*sr_it);

    transcode_node_to_global_dof_scheme(ncs_it, ncs_end, *sr_it);
  }

  entities_changed = true;
}

/* -------------------------------------------------------------------------- */
void DOFSynchronizer::fillEntityToSend(Array<UInt> & dofs_to_send) {
  UInt nb_dofs = dof_manager.getLocalSystemSize();

  this->entities_from_root.clear();
  dofs_to_send.resize(0);

  for (UInt d : arange(nb_dofs)) {
    if (not dof_manager.isLocalOrMasterDOF(d))
      continue;

    entities_from_root.push_back(d);
  }

  for (auto d : entities_from_root) {
    UInt global_dof = dof_manager.localToGlobalEquationNumber(d);
    dofs_to_send.push_back(global_dof);
  }
}

/* -------------------------------------------------------------------------- */
void DOFSynchronizer::onNodesAdded(const Array<UInt> & /*nodes_list*/) {
  auto dof_ids = dof_manager.getDOFIDs();

  for (auto sr : iterate_send_recv) {
    for (auto && data : communications.iterateSchemes(sr)) {
      auto & scheme = data.second;
      scheme.resize(0);
    }
  }

  for (auto & dof_id : dof_ids) {
    registerDOFs(dof_id);
  }

  // const auto & node_synchronizer =
  // dof_manager.getMesh().getNodeSynchronizer(); const auto &
  // node_communications = node_synchronizer.getCommunications();

  // std::map<UInt, std::vector<UInt>> nodes_per_proc[2];

  // for (auto sr : iterate_send_recv) {
  //   for (auto && data : node_communications.iterateSchemes(sr)) {
  //     auto proc = data.first;
  //     const auto & scheme = data.second;
  //     for (auto node : scheme) {
  //       nodes_per_proc[sr][proc].push_back(node);
  //     }
  //   }
  // }

  // std::map<UInt, std::vector<UInt>> dofs_per_proc[2];
  // for (auto & dof_id : dof_ids) {
  //   const auto & associated_nodes =
  //   dof_manager.getDOFsAssociatedNodes(dof_id); const auto &
  //   local_equation_numbers =
  //       dof_manager.getEquationsNumbers(dof_id);

  //   for (auto tuple : zip(associated_nodes, local_equation_numbers)) {
  //     UInt assoc_node;
  //     UInt local_eq_num;
  //     std::tie(assoc_node, local_eq_num) = tuple;

  //     for (auto sr_it = send_recv_t::begin(); sr_it != send_recv_t::end();
  //          ++sr_it) {
  //       for (auto & pair : nodes_per_proc[*sr_it]) {
  //         if (std::find(pair.second.end(), pair.second.end(), assoc_node) !=
  //             pair.second.end()) {
  //           dofs_per_proc[*sr_it][pair.first].push_back(local_eq_num);
  //         }
  //       }
  //     }
  //   }
  // }

  // for (auto sr_it = send_recv_t::begin(); sr_it != send_recv_t::end();
  //      ++sr_it) {
  //   for (auto & pair : dofs_per_proc[*sr_it]) {
  //     std::sort(pair.second.begin(), pair.second.end(),
  //               [this](UInt la, UInt lb) -> bool {
  //                 UInt ga = dof_manager.localToGlobalEquationNumber(la);
  //                 UInt gb = dof_manager.localToGlobalEquationNumber(lb);
  //                 return ga < gb;
  //               });

  //     auto & scheme = communications.getScheme(pair.first, *sr_it);
  //     scheme.resize(0);
  //     for (auto leq : pair.second) {
  //       scheme.push_back(leq);
  //     }
  //   }
  // }
  this->entities_changed = true;
}

} // namespace akantu
