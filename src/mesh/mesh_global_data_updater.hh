/**
 * @file   mesh_global_data_updater.hh
 *
 * @author Nicolas Richart
 *
 * @date creation  Sat Mar 03 2018
 *
 * @brief interface for the global data updater
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
#ifndef __AKANTU_MESH_GLOBAL_DATA_UPDATER_HH__
#define __AKANTU_MESH_GLOBAL_DATA_UPDATER_HH__

namespace akantu {

class NewNodesEvent;
class NewElementsEvent;

class MeshGlobalDataUpdater {
public:
  virtual ~MeshGlobalDataUpdater() = default;

  virtual std::tuple<UInt, UInt>
  updateData(NewNodesEvent & /*nodes_event*/,
             NewElementsEvent & /*elements_event*/) {
    return std::make_tuple(0, 0);
  }
};

} // namespace akantu

#endif /* __AKANTU_MESH_GLOBAL_DATA_UPDATER_HH__ */
