/**
 * @file   element_info_per_processor.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Fri Mar 11 14:56:42 2016
 *
 * @brief
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
#include "element_info_per_processor.hh"
#include "distributed_synchronizer.hh"
#include "static_communicator.hh"
/* -------------------------------------------------------------------------- */
#include <map>
#include <iostream>
#include <algorithm>
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
ElementInfoPerProc::ElementInfoPerProc(DistributedSynchronizer & synchronizer,
                                       StaticCommunicator & communicator,
                                       UInt message_cnt, UInt root, Mesh & mesh,
                                       ElementType type)
    : MeshAccessor(mesh), synchronizer(synchronizer), comm(communicator),
      rank(communicator.whoAmI()), nb_proc(communicator.getNbProc()),
      root(root), type(type), mesh(mesh),message_count(message_cnt) {}

/* -------------------------------------------------------------------------- */
void ElementInfoPerProc::fillCommunicationScheme(
    const Array<UInt> & partition) {
  AKANTU_DEBUG_IN();

  Element element;
  element.type = this->type;
  element.kind = Mesh::getKind(this->type);

  Array<UInt>::const_scalar_iterator part = partition.begin();
  for (UInt lel = 0; lel < nb_local_element; ++lel, ++part) {
    UInt nb_send = *part;
    element.element = lel;
    element.ghost_type = _not_ghost;
    for (UInt p = 0; p < nb_send; ++p, ++part) {
      UInt proc = *part;

      AKANTU_DEBUG(dblAccessory, "Must send : " << element << " to proc "
                                                << proc);
      (synchronizer.send_element[proc]).push_back(element);
    }
  }

  part = partition.begin();
  for (UInt gel = 0; gel < nb_ghost_element; ++gel, ++part) {
    UInt proc = *part;
    element.element = gel;
    element.ghost_type = _ghost;
    AKANTU_DEBUG(dblAccessory, "Must recv : " << element << " from proc "
                                              << proc);
    synchronizer.recv_element[proc].push_back(element);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

__END_AKANTU__
