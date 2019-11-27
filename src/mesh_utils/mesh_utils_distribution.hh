/**
 * @file   mesh_utils_distribution.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Sat Apr 01 2017
 *
 * @brief  Mesh utils to distribute a mesh
 *
 * @section LICENSE
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
#include "aka_common.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MESH_UTILS_DISTRIBUTION_HH__
#define __AKANTU_MESH_UTILS_DISTRIBUTION_HH__

namespace akantu {
class Mesh;
class MeshPartition;
} // namespace akantu

namespace akantu {
namespace MeshUtilsDistribution {
  /// Master call to distribute a mesh in a centralized manner (the UInt is just
  /// to avoid some shitty access from the slave...)
  void distributeMeshCentralized(Mesh & mesh, UInt,
                                 const MeshPartition & partition);
  /// Slave call to distribute a mesh in a centralized manner
  void distributeMeshCentralized(Mesh & mesh, UInt root);
} // namespace MeshUtilsDistribution

} // namespace akantu

#endif /* __AKANTU_MESH_UTILS_DISTRIBUTION_HH__ */
