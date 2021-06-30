/**
 * @file   mesh_partition_mesh_data.hh
 *
 * @author Dana Christen <dana.christen@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  mesh partitioning based on data provided in the mesh
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

#ifndef AKANTU_MESH_PARTITION_MESH_DATA_HH_
#define AKANTU_MESH_PARTITION_MESH_DATA_HH_

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh_partition.hh"

/* -------------------------------------------------------------------------- */

namespace akantu {

class MeshPartitionMeshData : public MeshPartition {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MeshPartitionMeshData(Mesh & mesh, UInt spatial_dimension,
                        const ID & id = "MeshPartitionerMeshData"
                        );

  MeshPartitionMeshData(Mesh & mesh,
                        const ElementTypeMapArray<UInt> & mapping,
                        UInt spatial_dimension,
                        const ID & id = "MeshPartitionerMeshData");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void partitionate(
      UInt nb_part,
      const std::function<Int(const Element &, const Element &)> &edge_load_func =
          [](auto && /*unused*/, auto && /*unused*/) { return 1; },
      const std::function<Int(const Element &)> &vertex_load_func =
          [](auto && /*unused*/) { return 1; }) override;

  void reorder() override;

  void setPartitionMapping(const ElementTypeMapArray<UInt> & mapping);

  void setPartitionMappingFromMeshData(const std::string & data_name);

private:
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  const ElementTypeMapArray<UInt> * partition_mapping;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

} // namespace akantu

#endif /* AKANTU_MESH_PARTITION_MESH_DATA_HH_ */
