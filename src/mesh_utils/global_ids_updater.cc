/**
 * @file   global_ids_updater.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Fri Apr 13 2012
 * @date last modification: Fri Oct 02 2015
 *
 * @brief  Functions of the GlobalIdsUpdater
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
#include "global_ids_updater.hh"
#include "element_synchronizer.hh"
#include "mesh_utils.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

UInt GlobalIdsUpdater::updateGlobalIDs(UInt old_nb_nodes) {
  UInt total_nb_new_nodes = this->updateGlobalIDsLocally(old_nb_nodes);

  this->synchronizeGlobalIDs();
  return total_nb_new_nodes;
}

UInt GlobalIdsUpdater::updateGlobalIDsLocally(UInt old_nb_nodes) {
  UInt total_nb_new_nodes =
      MeshUtils::updateLocalMasterGlobalConnectivity(mesh, old_nb_nodes);
  return total_nb_new_nodes;
}

void GlobalIdsUpdater::synchronizeGlobalIDs() {
  this->synchronizer.synchronizeOnce(*this, _gst_giu_global_conn);
}

} // akantu
