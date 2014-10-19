/**
 * @file   mesh_data.cc
 *
 * @author Dana Christen <dana.christen@gmail.com>
 *
 * @date creation: Fri May 03 2013
 * @date last modification: Thu Jun 05 2014
 *
 * @brief  Stores generic data loaded from the mesh file
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#include "mesh.hh"
#include "mesh_data.hh"

__BEGIN_AKANTU__

MeshData::MeshData(const ID & _id, const ID & parent_id, const MemoryID & mem_id)
  : Memory(parent_id + ":" + _id, mem_id) {
}

MeshData::~MeshData() {
  ElementalDataMap::iterator it, end;
  for(it = elemental_data.begin(); it != elemental_data.end(); ++it) {
    delete it->second;
  }
}

__END_AKANTU__

