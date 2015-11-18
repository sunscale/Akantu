/**
 * @file   global_ids_updater.cc
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @date   Fri Oct  2 13:44:02 2015
 *
 * @brief  Functions of the GlobalIdsUpdater
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
#include "global_ids_updater.hh"
#include "mesh_utils.hh"

__BEGIN_AKANTU__

UInt GlobalIdsUpdater::updateGlobalIDs(UInt old_nb_nodes) {
  UInt total_nb_new_nodes = MeshUtils::updateLocalMasterGlobalConnectivity(mesh, old_nb_nodes);

  synchronizer->computeBufferSize(*this, _gst_giu_global_conn);
  synchronizer->asynchronousSynchronize(*this, _gst_giu_global_conn);
  synchronizer->waitEndSynchronize(*this, _gst_giu_global_conn);

  return total_nb_new_nodes;
}

__END_AKANTU__