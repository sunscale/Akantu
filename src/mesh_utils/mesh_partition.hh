/**
 * @file   mesh_partition.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Tue Jan 23 2018
 *
 * @brief  tools to partitionate a mesh
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
#include "aka_csr.hh"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MESH_PARTITION_HH_
#define AKANTU_MESH_PARTITION_HH_


namespace akantu {

class MeshPartition {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MeshPartition(Mesh & mesh, UInt spatial_dimension,
                const ID & id = "MeshPartitioner");

  virtual ~MeshPartition();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// define a partition of the mesh
  virtual void partitionate(
      UInt nb_part,
      const std::function<Int(const Element &, const Element &)> &
          edge_load_func =
              [](auto && /*unused*/, auto && /*unused*/) { return 1; },
      const std::function<Int(const Element &)> & vertex_load_func =
          [](auto && /*unused*/) { return 1; }) = 0;

  /// reorder the nodes to reduce the filling during the factorization of a
  /// matrix that has a profil based on the connectivity of the mesh
  virtual void reorder() = 0;

  /// fill the partitions array with a given linearized partition information
  void fillPartitionInformation(const Mesh & mesh,
                                const Int * linearized_partitions);

  virtual void printself(std::ostream & stream, int indent = 0) const;

protected:
  /// build the dual graph of the mesh, for all element of spatial_dimension
  void
  buildDualGraph(Array<Int> & dxadj, Array<Int> & dadjncy,
                 Array<Int> & edge_loads,
                 const std::function<Int(const Element &, const Element &)> &
                     edge_load_func,
                 Array<Int> & vertex_loads,
                 const std::function<Int(const Element &)> & vertex_load_func);

  /// tweak the mesh to handle the PBC pairs
  void tweakConnectivity();
  /// restore the mesh that has been tweaked
  void restoreConnectivity();

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  bool hasPartitions(ElementType type, GhostType ghost_type);
  AKANTU_GET_MACRO(Partitions, partitions, const ElementTypeMapArray<UInt> &);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Partition, partitions, UInt);

  AKANTU_GET_MACRO(GhostPartitionCSR, ghost_partitions_csr,
                   const ElementTypeMap<CSR<UInt>> &);

  AKANTU_GET_MACRO(NbPartition, nb_partitions, UInt);
  AKANTU_SET_MACRO(NbPartition, nb_partitions, UInt);

protected:
  UInt linearized(const Element & element);
  Element unlinearized(UInt lin_element);
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// id
  ID id;

  /// the mesh to partition
  Mesh & mesh;

  /// dimension of the elements to consider in the mesh
  UInt spatial_dimension;

  /// number of partitions
  UInt nb_partitions;

  /// partition numbers
  ElementTypeMapArray<UInt> partitions;

  ElementTypeMap<CSR<UInt>> ghost_partitions_csr;
  ElementTypeMapArray<UInt> ghost_partitions;
  ElementTypeMapArray<UInt> ghost_partitions_offset;

  Array<UInt> * permutation;

  ElementTypeMapArray<UInt> saved_connectivity;

  // vector of pair to ensure the iteration order
  std::vector<std::pair<ElementType, UInt>> linearized_offsets;
};

/// standard output stream operator
inline std::ostream & operator<<(std::ostream & stream,
                                 const MeshPartition & _this) {
  _this.printself(stream);
  return stream;
}

} // namespace akantu

#ifdef AKANTU_USE_SCOTCH
#include "mesh_partition_scotch.hh"
#endif

#endif /* AKANTU_MESH_PARTITION_HH_ */
