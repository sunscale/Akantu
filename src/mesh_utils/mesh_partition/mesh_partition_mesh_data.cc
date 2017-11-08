/**
 * @file   mesh_partition_mesh_data.cc
 *
 * @author Dana Christen <dana.christen@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Fri May 03 2013
 * @date last modification: Wed Nov 11 2015
 *
 * @brief  implementation of the MeshPartitionMeshData class
 *
 * @section LICENSE
 *
 * Copyright  (©)  2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de Lausanne)
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

/* -------------------------------------------------------------------------- */
#include "mesh_partition_mesh_data.hh"
#if !defined(AKANTU_NDEBUG)
#include <set>
#endif

/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
MeshPartitionMeshData::MeshPartitionMeshData(const Mesh & mesh,
                                             UInt spatial_dimension,
                                             const ID & id,
                                             const MemoryID & memory_id)
    : MeshPartition(mesh, spatial_dimension, id, memory_id) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
MeshPartitionMeshData::MeshPartitionMeshData(
    const Mesh & mesh, const ElementTypeMapArray<UInt> & mapping,
    UInt spatial_dimension, const ID & id, const MemoryID & memory_id)
    : MeshPartition(mesh, spatial_dimension, id, memory_id),
      partition_mapping(&mapping) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshPartitionMeshData::partitionate(UInt nb_part,
                                         __attribute__((unused))
                                         const EdgeLoadFunctor & edge_load_func,
                                         const Array<UInt> & pairs) {
  AKANTU_DEBUG_IN();

  tweakConnectivity(pairs);

  nb_partitions = nb_part;
  GhostType ghost_type = _not_ghost;
  UInt spatial_dimension = mesh.getSpatialDimension();
  Mesh::type_iterator it =
      mesh.firstType(spatial_dimension, ghost_type, _ek_not_defined);
  Mesh::type_iterator end =
      mesh.lastType(spatial_dimension, ghost_type, _ek_not_defined);

  UInt linearized_el = 0;
  UInt nb_elements = mesh.getNbElement(mesh.getSpatialDimension(), ghost_type);
  auto * partition_list = new Int[nb_elements];

#if !defined(AKANTU_NDEBUG)
  std::set<UInt> partitions;
#endif
  for (; it != end; ++it) {
    ElementType type = *it;
    const Array<UInt> & partition_array =
        (*partition_mapping)(type, ghost_type);
    Array<UInt>::const_iterator<Vector<UInt> > p_it = partition_array.begin(1);
    Array<UInt>::const_iterator<Vector<UInt> > p_end = partition_array.end(1);
    AKANTU_DEBUG_ASSERT(UInt(p_end - p_it) ==
                            mesh.getNbElement(type, ghost_type),
                        "The partition mapping does not have the right number "
                            << "of entries for type " << type
                            << " and ghost type " << ghost_type << "."
                            << " Tags=" << p_end - p_it
                            << " Mesh=" << mesh.getNbElement(type, ghost_type));
    for (; p_it != p_end; ++p_it, ++linearized_el) {
      partition_list[linearized_el] = (*p_it)(0);
#if !defined(AKANTU_NDEBUG)
      partitions.insert((*p_it)(0));
#endif
    }
  }

#if !defined(AKANTU_NDEBUG)
  AKANTU_DEBUG_ASSERT(partitions.size() == nb_part,
                      "The number of real partitions does not match with the "
                      "number of asked partitions");
#endif

  fillPartitionInformation(mesh, partition_list);

  delete[] partition_list;

  restoreConnectivity();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshPartitionMeshData::reorder() { AKANTU_DEBUG_TO_IMPLEMENT(); }

/* -------------------------------------------------------------------------- */
void MeshPartitionMeshData::setPartitionMapping(
    const ElementTypeMapArray<UInt> & mapping) {
  partition_mapping = &mapping;
}

/* -------------------------------------------------------------------------- */
void MeshPartitionMeshData::setPartitionMappingFromMeshData(
    const std::string & data_name) {
  partition_mapping = &(mesh.getData<UInt>(data_name));
}

} // akantu
