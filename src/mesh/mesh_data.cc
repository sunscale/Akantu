/**
 * @file   mesh_data.cc
 *
 * @author Dana Christen <dana.christen@gmail.com>
 *
 * @date creation: Fri Apr 13 2012
 * @date last modification: Mon Jun 19 2017
 *
 * @brief  Stores generic data loaded from the mesh file
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

#include "mesh_data.hh"
#include "mesh.hh"

namespace akantu {

MeshData::MeshData(const ID & _id, const ID & parent_id)
    : _id(parent_id + ":" + _id) {}

} // namespace akantu
