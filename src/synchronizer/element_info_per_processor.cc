/**
 * @file   element_info_per_processor.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Mar 16 2016
 * @date last modification: Tue Nov 07 2017
 *
 * @brief  Helper class to distribute a mesh
 *
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
#include "element_info_per_processor.hh"
#include "communicator.hh"
#include "element_synchronizer.hh"
/* -------------------------------------------------------------------------- */
#include <algorithm>
#include <iostream>
#include <map>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
ElementInfoPerProc::ElementInfoPerProc(ElementSynchronizer & synchronizer,
                                       UInt message_cnt, UInt root,
                                       ElementType type)
    : MeshAccessor(synchronizer.getMesh()), synchronizer(synchronizer),
      rank(synchronizer.getCommunicator().whoAmI()),
      nb_proc(synchronizer.getCommunicator().getNbProc()), root(root),
      type(type), message_count(message_cnt), mesh(synchronizer.getMesh()),
      comm(synchronizer.getCommunicator()) {}

/* -------------------------------------------------------------------------- */
bool ElementInfoPerProc::synchronize() {
  auto need_synchronize = needSynchronize();

  if (need_synchronize) {
    synchronizeConnectivities();
    synchronizePartitions();
    synchronizeTags();
    synchronizeGroups();
  }

  return need_synchronize;
}

/* -------------------------------------------------------------------------- */
void ElementInfoPerProc::fillCommunicationScheme(
    const Array<UInt> & partition) {
  AKANTU_DEBUG_IN();

  Element element;
  element.type = this->type;

  auto & communications = this->synchronizer.getCommunications();
  auto part = partition.begin();

  std::map<UInt, Array<Element>> send_array_per_proc;
  for (UInt lel = 0; lel < nb_local_element; ++lel) {
    UInt nb_send = *part;
    ++part;

    element.element = lel;
    element.ghost_type = _not_ghost;
    for (UInt p = 0; p < nb_send; ++p, ++part) {
      UInt proc = *part;

      AKANTU_DEBUG(dblAccessory,
                   "Must send : " << element << " to proc " << proc);
      send_array_per_proc[proc].push_back(element);
    }
  }

  for (auto & send_schemes : send_array_per_proc) {
    if (send_schemes.second.size() == 0)
      continue;
    auto & scheme = communications.createSendScheme(send_schemes.first);
    scheme.append(send_schemes.second);
  }

  std::map<UInt, Array<Element>> recv_array_per_proc;

  for (UInt gel = 0; gel < nb_ghost_element; ++gel, ++part) {
    UInt proc = *part;
    element.element = gel;
    element.ghost_type = _ghost;
    AKANTU_DEBUG(dblAccessory,
                 "Must recv : " << element << " from proc " << proc);
    recv_array_per_proc[proc].push_back(element);
  }

  for (auto & recv_schemes : recv_array_per_proc) {
    if (recv_schemes.second.size() == 0)
      continue;
    auto & scheme = communications.createRecvScheme(recv_schemes.first);
    scheme.append(recv_schemes.second);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

} // namespace akantu
