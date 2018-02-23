/**
 * @file   mesh_partition.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Aug 17 2010
 * @date last modification: Fri Jan 22 2016
 *
 * @brief  implementation of common part of all partitioner
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
#include "mesh_partition.hh"
#include "aka_iterators.hh"
#include "aka_types.hh"
#include "mesh_accessor.hh"
#include "mesh_utils.hh"
/* -------------------------------------------------------------------------- */
#include <algorithm>
#include <numeric>
#include <unordered_map>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
MeshPartition::MeshPartition(const Mesh & mesh, UInt spatial_dimension,
                             const ID & id, const MemoryID & memory_id)
    : Memory(id, memory_id), mesh(mesh), spatial_dimension(spatial_dimension),
      partitions("partition", id, memory_id),
      ghost_partitions("ghost_partition", id, memory_id),
      ghost_partitions_offset("ghost_partition_offset", id, memory_id),
      saved_connectivity("saved_connectivity", id, memory_id) {
  AKANTU_DEBUG_IN();

  UInt nb_total_element = 0;
  for (auto && type :
       mesh.elementTypes(spatial_dimension, _not_ghost, _ek_not_defined)) {
    linearized_offsets.push_back(std::make_pair(type, nb_total_element));
    nb_total_element += mesh.getConnectivity(type).size();
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
MeshPartition::~MeshPartition() = default;

/* -------------------------------------------------------------------------- */
UInt MeshPartition::linearized(const Element & element) {
  auto it =
      std::find_if(linearized_offsets.begin(), linearized_offsets.end(),
                   [&element](auto & a) { return a.first == element.type; });
  AKANTU_DEBUG_ASSERT(it != linearized_offsets.end(),
                      "A bug might be crawling around this corner...");
  return (it->second + element.element);
}

/* -------------------------------------------------------------------------- */
Element MeshPartition::unlinearized(UInt lin_element) {
  ElementType type{_not_defined};
  UInt offset{0};
  for (auto & pair : linearized_offsets) {
    if (lin_element < pair.second)
      continue;
    std::tie(type, offset) = pair;
  }

  return Element{type, lin_element - offset, _not_ghost};
}

/* -------------------------------------------------------------------------- */
/**
 * TODO this function should most probably be rewritten in a more modern way
 * conversion in c++ of the GENDUALMETIS (mesh.c) function wrote by George in
 * Metis (University of Minnesota)
 */
void MeshPartition::buildDualGraph(Array<Int> & dxadj, Array<Int> & dadjncy,
                                   Array<Int> & edge_loads,
                                   const EdgeLoadFunctor & edge_load_func) {
  AKANTU_DEBUG_IN();

  std::map<ElementType, std::tuple<const Array<UInt> *, UInt, UInt>>
      connectivities;
  UInt spatial_dimension = mesh.getSpatialDimension();
  UInt nb_total_element{0};

  for (auto & type :
       mesh.elementTypes(spatial_dimension, _not_ghost, _ek_not_defined)) {
    auto type_p1 = mesh.getP1ElementType(type);
    auto nb_nodes_per_element_p1 = mesh.getNbNodesPerElement(type_p1);
    auto magic_number = mesh.getNbNodesPerElement(mesh.getFacetType(type_p1));

    const auto & conn = mesh.getConnectivity(type, _not_ghost);
    connectivities[type] = std::make_tuple(
        &conn, nb_nodes_per_element_p1, magic_number);
    nb_total_element += conn.size();
  }

  CSR<Element> node_to_elem;
  MeshUtils::buildNode2Elements(mesh, node_to_elem);

  dxadj.resize(nb_total_element + 1);
  /// initialize the dxadj array
  auto dxadj_it = dxadj.begin();
  for (auto & pair : connectivities) {
    const auto & connectivity = *std::get<0>(pair.second);
    auto nb_nodes_per_element_p1 = std::get<1>(pair.second);

    std::fill_n(dxadj_it, connectivity.size(), nb_nodes_per_element_p1);
    dxadj_it += connectivity.size();
  }

  /// convert the dxadj_val array in a csr one
  for (UInt i = 1; i < nb_total_element; ++i)
    dxadj(i) += dxadj(i - 1);
  for (UInt i = nb_total_element; i > 0; --i)
    dxadj(i) = dxadj(i - 1);
  dxadj(0) = 0;

  dadjncy.resize(2 * dxadj(nb_total_element));
  edge_loads.resize(2 * dxadj(nb_total_element));

  /// weight map to determine adjacency
  std::unordered_map<UInt, UInt> weight_map;

  for (auto & pair : connectivities) {
    auto type = pair.first;
    const auto & connectivity = *std::get<0>(pair.second);
    auto nb_nodes_per_element = std::get<1>(pair.second);
    auto magic_number = std::get<2>(pair.second);

    Element element{type, 0, _not_ghost};

    for (const auto & conn :
         make_view(connectivity, connectivity.getNbComponent())) {
      auto linearized_el = linearized(element);

      /// fill the weight map
      for (UInt n : arange(nb_nodes_per_element)) {
        auto && node = conn(n);
        for (auto k = node_to_elem.rbegin(node); k != node_to_elem.rend(node);
             --k) {
          auto & current_element = *k;
          auto current_el = linearized(current_element);
          AKANTU_DEBUG_ASSERT(current_el != UInt(-1),
                              "Linearized element not found");
          if (current_el <= linearized_el)
            break;

          auto weight_map_insert =
              weight_map.insert(std::make_pair(current_el, 1));
          if (not weight_map_insert.second)
            (weight_map_insert.first->second)++;
        }
      }

      /// each element with a weight of the size of a facet are adjacent
      for (auto & weight_pair : weight_map) {
        auto & adjacent_el = weight_pair.first;
        auto magic = weight_pair.second;
        if (magic != magic_number)
          continue;

#if defined(AKANTU_COHESIVE_ELEMENT)
        /// Patch in order to prevent neighboring cohesive elements
        /// from detecting each other
        auto adjacent_element = unlinearized(adjacent_el);

        auto el_kind = element.kind();
        auto adjacent_el_kind = adjacent_element.kind();

        if (el_kind == adjacent_el_kind && el_kind == _ek_cohesive)
          continue;
#endif
        UInt index_adj = dxadj(adjacent_el)++;
        UInt index_lin = dxadj(linearized_el)++;

        dadjncy(index_lin) = adjacent_el;
        dadjncy(index_adj) = linearized_el;
      }

      element.element++;
      weight_map.clear();
    }
  }

  Int k_start = 0, linerized_el = 0, j = 0;
  for (auto & pair : connectivities) {
    const auto & connectivity = *std::get<0>(pair.second);
    auto nb_nodes_per_element_p1 = std::get<1>(pair.second);
    auto nb_element = connectivity.size();

    for (UInt el = 0; el < nb_element; ++el, ++linerized_el) {
      for (Int k = k_start; k < dxadj(linerized_el); ++k, ++j)
        dadjncy(j) = dadjncy(k);
      dxadj(linerized_el) = j;
      k_start += nb_nodes_per_element_p1;
    }
  }

  for (UInt i = nb_total_element; i > 0; --i)
    dxadj(i) = dxadj(i - 1);
  dxadj(0) = 0;

  UInt adj = 0;
  for (UInt i = 0; i < nb_total_element; ++i) {
    UInt nb_adj = dxadj(i + 1) - dxadj(i);
    for (UInt j = 0; j < nb_adj; ++j, ++adj) {
      Int el_adj_id = dadjncy(dxadj(i) + j);
      Element el = unlinearized(i);
      Element el_adj = unlinearized(el_adj_id);

      Int load = edge_load_func(el, el_adj);
      edge_loads(adj) = load;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshPartition::fillPartitionInformation(
    const Mesh & mesh, const Int * linearized_partitions) {
  AKANTU_DEBUG_IN();

  CSR<Element> node_to_elem;

  MeshUtils::buildNode2Elements(mesh, node_to_elem);

  UInt linearized_el = 0;
  for (auto & type :
       mesh.elementTypes(spatial_dimension, _not_ghost, _ek_not_defined)) {
    UInt nb_element = mesh.getNbElement(type);
    UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

    partitions.alloc(nb_element, 1, type, _not_ghost);
    auto & ghost_part_csr = ghost_partitions_csr(type, _not_ghost);
    ghost_part_csr.resizeRows(nb_element);

    ghost_partitions_offset.alloc(nb_element + 1, 1, type, _ghost);
    ghost_partitions.alloc(0, 1, type, _ghost);

    const Array<UInt> & connectivity = mesh.getConnectivity(type, _not_ghost);

    for (UInt el = 0; el < nb_element; ++el, ++linearized_el) {
      UInt part = linearized_partitions[linearized_el];

      partitions(type, _not_ghost)(el) = part;
      std::list<UInt> list_adj_part;
      for (UInt n = 0; n < nb_nodes_per_element; ++n) {
        UInt node = connectivity.storage()[el * nb_nodes_per_element + n];
        CSR<Element>::iterator ne;
        for (ne = node_to_elem.begin(node); ne != node_to_elem.end(node);
             ++ne) {
          const Element & adj_element = *ne;
          UInt adj_el = linearized(adj_element);
          UInt adj_part = linearized_partitions[adj_el];
          if (part != adj_part) {
            list_adj_part.push_back(adj_part);
          }
        }
      }

      list_adj_part.sort();
      list_adj_part.unique();

      for (auto & adj_part : list_adj_part) {
        ghost_part_csr.getRows().push_back(adj_part);
        ghost_part_csr.rowOffset(el)++;

        ghost_partitions(type, _ghost).push_back(adj_part);
        ghost_partitions_offset(type, _ghost)(el)++;
      }
    }

    ghost_part_csr.countToCSR();

    /// convert the ghost_partitions_offset array in an offset array
    Array<UInt> & ghost_partitions_offset_ptr =
        ghost_partitions_offset(type, _ghost);
    for (UInt i = 1; i < nb_element; ++i)
      ghost_partitions_offset_ptr(i) += ghost_partitions_offset_ptr(i - 1);
    for (UInt i = nb_element; i > 0; --i)
      ghost_partitions_offset_ptr(i) = ghost_partitions_offset_ptr(i - 1);
    ghost_partitions_offset_ptr(0) = 0;
  }

  // All Facets
  for (Int sp = spatial_dimension - 1; sp >= 0; --sp) {
    for (auto & type : mesh.elementTypes(sp, _not_ghost, _ek_not_defined)) {
      UInt nb_element = mesh.getNbElement(type);

      partitions.alloc(nb_element, 1, type, _not_ghost);
      AKANTU_DEBUG_INFO("Allocating partitions for " << type);
      auto & ghost_part_csr = ghost_partitions_csr(type, _not_ghost);
      ghost_part_csr.resizeRows(nb_element);

      ghost_partitions_offset.alloc(nb_element + 1, 1, type, _ghost);
      ghost_partitions.alloc(0, 1, type, _ghost);
      AKANTU_DEBUG_INFO("Allocating ghost_partitions for " << type);
      const Array<std::vector<Element>> & elem_to_subelem =
          mesh.getElementToSubelement(type, _not_ghost);
      for (UInt i(0); i < mesh.getNbElement(type, _not_ghost);
           ++i) { // Facet loop

        const std::vector<Element> & adjacent_elems = elem_to_subelem(i);
        if (!adjacent_elems.empty()) {
          Element min_elem{_max_element_type, std::numeric_limits<UInt>::max(),
                           *ghost_type_t::end()};
          UInt min_part(std::numeric_limits<UInt>::max());
          std::set<UInt> adjacent_parts;

          for (UInt j(0); j < adjacent_elems.size(); ++j) {
            UInt adjacent_elem_id = adjacent_elems[j].element;
            UInt adjacent_elem_part =
                partitions(adjacent_elems[j].type,
                           adjacent_elems[j].ghost_type)(adjacent_elem_id);
            if (adjacent_elem_part < min_part) {
              min_part = adjacent_elem_part;
              min_elem = adjacent_elems[j];
            }
            adjacent_parts.insert(adjacent_elem_part);
          }
          partitions(type, _not_ghost)(i) = min_part;

          auto git = ghost_partitions_csr(min_elem.type, _not_ghost)
                         .begin(min_elem.element);
          auto gend = ghost_partitions_csr(min_elem.type, _not_ghost)
                          .end(min_elem.element);
          for (; git != gend; ++git) {

            adjacent_parts.insert(min_part);
          }
          adjacent_parts.erase(min_part);
          for (auto & part : adjacent_parts) {
            ghost_part_csr.getRows().push_back(part);
            ghost_part_csr.rowOffset(i)++;
            ghost_partitions(type, _ghost).push_back(part);
          }

          ghost_partitions_offset(type, _ghost)(i + 1) =
              ghost_partitions_offset(type, _ghost)(i + 1) +
              adjacent_elems.size();
        } else {
          partitions(type, _not_ghost)(i) = 0;
        }
      }
      ghost_part_csr.countToCSR();
    }
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshPartition::tweakConnectivity(const Array<UInt> & pairs) {
  AKANTU_DEBUG_IN();

  if (pairs.size() == 0)
    return;

  Mesh::type_iterator it =
      mesh.firstType(spatial_dimension, _not_ghost, _ek_not_defined);
  Mesh::type_iterator end =
      mesh.lastType(spatial_dimension, _not_ghost, _ek_not_defined);

  for (; it != end; ++it) {
    ElementType type = *it;

    Array<UInt> & conn =
        const_cast<Array<UInt> &>(mesh.getConnectivity(type, _not_ghost));
    UInt nb_nodes_per_element = conn.getNbComponent();
    UInt nb_element = conn.size();

    Array<UInt> & saved_conn = saved_connectivity.alloc(
        nb_element, nb_nodes_per_element, type, _not_ghost);
    saved_conn.copy(conn);

    for (UInt i = 0; i < pairs.size(); ++i) {
      for (UInt el = 0; el < nb_element; ++el) {
        for (UInt n = 0; n < nb_nodes_per_element; ++n) {
          if (pairs(i, 1) == conn(el, n))
            conn(el, n) = pairs(i, 0);
        }
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshPartition::restoreConnectivity() {
  AKANTU_DEBUG_IN();
  MeshAccessor mesh_accessor(const_cast<Mesh &>(mesh));
  for (auto && type : saved_connectivity.elementTypes(
           spatial_dimension, _not_ghost, _ek_not_defined)) {
    auto & conn = mesh_accessor.getConnectivity(type, _not_ghost);
    auto & saved_conn = saved_connectivity(type, _not_ghost);
    conn.copy(saved_conn);
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
bool MeshPartition::hasPartitions(const ElementType & type,
                                  const GhostType & ghost_type) {
  return partitions.exists(type, ghost_type);
}

/* -------------------------------------------------------------------------- */

} // namespace akantu
