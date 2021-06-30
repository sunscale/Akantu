/**
 * @file   mesh_partition_scotch.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  mesh partitioning based on libScotch
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
#include "mesh_partition.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MESH_PARTITION_SCOTCH_HH_
#define AKANTU_MESH_PARTITION_SCOTCH_HH_

/* -------------------------------------------------------------------------- */

namespace akantu {

class MeshPartitionScotch : public MeshPartition {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MeshPartitionScotch(Mesh & mesh, UInt spatial_dimension,
                      const ID & id = "mesh_partition_scotch");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void partitionate(
      UInt nb_part,
      const std::function<Int(const Element &, const Element &)> &
          edge_load_func =
              [](auto && /*unused*/, auto && /*unused*/) { return 1; },
      const std::function<Int(const Element &)> & vertex_load_func =
          [](auto && /*unused*/) { return 1; }) override;

  void reorder() override;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
};

} // namespace akantu

#endif /* AKANTU_MESH_PARTITION_SCOTCH_HH_ */
