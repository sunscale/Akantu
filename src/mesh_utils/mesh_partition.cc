/**
 * @file   mesh_partition.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Aug 17 2010
 * @date last modification: Wed Jan 24 2018
 *
 * @brief  implementation of common part of all partitioner
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
#include "aka_iterators.hh"
#include "aka_types.hh"
#include "mesh_accessor.hh"
#include "mesh_iterators.hh"
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
 * conversion in c++ of the METIS_MeshToDual (mesh.c) function wrote by George
 * in Metis (University of Minnesota)
 */
void MeshPartition::buildDualGraph(
    Array<Int> & dxadj, Array<Int> & dadjncy, Array<Int> & edge_loads,
    std::function<Int(const Element &, const Element &)> edge_load_func,
    Array<Int> & vertex_loads,
    std::function<Int(const Element &)> vertex_load_func) {

  CSR<Element> nodes_to_elements;
  MeshUtils::buildNode2Elements(mesh, nodes_to_elements);

  std::unordered_map<UInt, std::vector<UInt>> adjacent_elements;

  // for each elements look for its connected elements
  for_each_element(
      mesh,
      [&](auto && element) {
        const auto & conn = mesh.getConnectivity(element);
        std::map<Element, UInt> hits;

        // count the number of nodes shared with a given element
        for (auto && node : conn) {
          for (auto && connected_element : nodes_to_elements.getRow(node)) {
            ++hits[connected_element];
          }
        }

        // define a minumum number of nodes to share to be considered as a
        // ajacent element
        UInt magic_number{conn.size()};
        for (auto n : arange(mesh.getNbFacetTypes(element.type))) {
          magic_number = std::min(
              mesh.getNbNodesPerElement(mesh.getFacetType(element.type, n)),
              magic_number);
        }

        // check all neighbors to see which ones are "adjacent"
        for (auto && data : hits) {
          const auto & adjacent_element = data.first;
          // not adjacent to miself
          if (adjacent_element == element)
            continue;

          // not enough shared nodes
          if (data.second < magic_number)
            continue;

        /// Patch in order to prevent neighboring cohesive elements
        /// from detecting each other
#if defined(AKANTU_COHESIVE_ELEMENT)
          auto element_kind = element.kind();
          auto adjacent_element_kind = adjacent_element.kind();

          if (element_kind == adjacent_element_kind &&
              element_kind == _ek_cohesive)
            continue;
#endif

          adjacent_elements[this->linearized(element)].push_back(
              this->linearized(adjacent_element));
        }
      },
      _spatial_dimension = mesh.getSpatialDimension(),
      _element_kind = _ek_not_defined);

  // prepare the arrays
  auto nb_elements{adjacent_elements.size()};
  dxadj.resize(nb_elements + 1);
  vertex_loads.resize(nb_elements);

  for (auto && data : adjacent_elements) {
    const auto & element{data.first};
    const auto & neighbors{data.second};
    dxadj[element] = neighbors.size();
  }

  /// convert the dxadj array of sizes in a csr one of offsets
  for (UInt i = 1; i < nb_elements; ++i)
    dxadj(i) += dxadj(i - 1);
  for (UInt i = nb_elements; i > 0; --i)
    dxadj(i) = dxadj(i - 1);
  dxadj(0) = 0;

  dadjncy.resize(dxadj(nb_elements));
  edge_loads.resize(dadjncy.size());

  // fill the different arrays
  for (auto && data : adjacent_elements) {
    const auto & element{data.first};
    const auto & neighbors{data.second};

    auto unlinearized_element = unlinearized(element);
    vertex_loads(element) = vertex_load_func(unlinearized_element);

    auto pos = dxadj(element);

    for (auto && neighbor : neighbors) {
      dadjncy(pos) = neighbor;
      edge_loads(pos) =
          edge_load_func(unlinearized_element, unlinearized(neighbor));
      ++pos;
    }
  }
}

/* -------------------------------------------------------------------------- */
/**
 * TODO this function should most probably be rewritten in a more modern way
 * conversion in c++ of the GENDUALMETIS (mesh.c) function wrote by George in
 * Metis (University of Minnesota)
 */
// void MeshPartition::buildDualGraph(
//     Array<Int> & dxadj, Array<Int> & dadjncy, Array<Int> & edge_loads,
//     std::function<Int(const Element &, const Element &)> edge_load_func,
//     Array<Int> & vertex_loads,
//     std::function<Int(const Element &)> vertex_load_func) {
//   AKANTU_DEBUG_IN();

//   std::map<ElementType, std::tuple<const Array<UInt> *, UInt, UInt>>
//       connectivities;
//   UInt spatial_dimension = mesh.getSpatialDimension();
//   UInt nb_total_element{0};

//   for (auto & type :
//        mesh.elementTypes(spatial_dimension, _not_ghost, _ek_not_defined)) {
//     auto type_p1 = mesh.getP1ElementType(type);
//     auto nb_nodes_per_element_p1 = mesh.getNbNodesPerElement(type_p1);
//     const auto & conn = mesh.getConnectivity(type, _not_ghost);

//     for (auto n : arange(mesh.getNbFacetTypes(type_p1))) {
//       auto magic_number =
//           mesh.getNbNodesPerElement(mesh.getFacetType(type_p1, n));
//       connectivities[type] =
//           std::make_tuple(&conn, nb_nodes_per_element_p1, magic_number);
//     }

//     nb_total_element += conn.size();
//   }

//   CSR<Element> node_to_elem;
//   MeshUtils::buildNode2Elements(mesh, node_to_elem);

//   dxadj.resize(nb_total_element + 1);
//   /// initialize the dxadj array
//   auto dxadj_it = dxadj.begin();
//   for (auto & pair : connectivities) {
//     const auto & connectivity = *std::get<0>(pair.second);
//     auto nb_nodes_per_element_p1 = std::get<1>(pair.second);

//     std::fill_n(dxadj_it, connectivity.size(), nb_nodes_per_element_p1);
//     dxadj_it += connectivity.size();
//   }

//   /// convert the dxadj_val array in a csr one
//   for (UInt i = 1; i < nb_total_element; ++i)
//     dxadj(i) += dxadj(i - 1);
//   for (UInt i = nb_total_element; i > 0; --i)
//     dxadj(i) = dxadj(i - 1);
//   dxadj(0) = 0;

//   dadjncy.resize(2 * dxadj(nb_total_element));

//   /// weight map to determine adjacency
//   std::unordered_map<UInt, UInt> weight_map;

//   for (auto & pair : connectivities) {
//     auto type = pair.first;
//     const auto & connectivity = *std::get<0>(pair.second);
//     auto nb_nodes_per_element = std::get<1>(pair.second);
//     auto magic_number = std::get<2>(pair.second);

//     Element element{type, 0, _not_ghost};

//     for (const auto & conn :
//          make_view(connectivity, connectivity.getNbComponent())) {
//       auto linearized_el = linearized(element);

//       /// fill the weight map
//       for (UInt n : arange(nb_nodes_per_element)) {
//         auto && node = conn(n);

//         for (auto k = node_to_elem.rbegin(node); k !=
//         node_to_elem.rend(node);
//              --k) {
//           auto & current_element = *k;
//           auto current_el = linearized(current_element);
//           AKANTU_DEBUG_ASSERT(current_el != UInt(-1),
//                               "Linearized element not found");
//           if (current_el <= linearized_el)
//             break;

//           auto weight_map_insert =
//               weight_map.insert(std::make_pair(current_el, 1));
//           if (not weight_map_insert.second)
//             (weight_map_insert.first->second)++;
//         }
//       }

//       /// each element with a weight of the size of a facet are adjacent
//       for (auto & weight_pair : weight_map) {
//         auto & adjacent_el = weight_pair.first;
//         auto magic = weight_pair.second;
//         if (magic != magic_number)
//           continue;

// #if defined(AKANTU_COHESIVE_ELEMENT)
//         /// Patch in order to prevent neighboring cohesive elements
//         /// from detecting each other
//         auto adjacent_element = unlinearized(adjacent_el);

//         auto el_kind = element.kind();
//         auto adjacent_el_kind = adjacent_element.kind();

//         if (el_kind == adjacent_el_kind && el_kind == _ek_cohesive)
//           continue;
// #endif
//         UInt index_adj = dxadj(adjacent_el)++;
//         UInt index_lin = dxadj(linearized_el)++;

//         dadjncy(index_lin) = adjacent_el;
//         dadjncy(index_adj) = linearized_el;
//       }

//       element.element++;
//       weight_map.clear();
//     }
//   }

//   Int k_start = 0, linerized_el = 0, j = 0;
//   for (auto & pair : connectivities) {
//     const auto & connectivity = *std::get<0>(pair.second);
//     auto nb_nodes_per_element_p1 = std::get<1>(pair.second);
//     auto nb_element = connectivity.size();

//     for (UInt el = 0; el < nb_element; ++el, ++linerized_el) {
//       for (Int k = k_start; k < dxadj(linerized_el); ++k, ++j)
//         dadjncy(j) = dadjncy(k);
//       dxadj(linerized_el) = j;
//       k_start += nb_nodes_per_element_p1;
//     }
//   }

//   for (UInt i = nb_total_element; i > 0; --i)
//     dxadj(i) = dxadj(i - 1);
//   dxadj(0) = 0;

//   vertex_loads.resize(dxadj.size() - 1);
//   edge_loads.resize(dadjncy.size());
//   UInt adj = 0;
//   for (UInt i = 0; i < nb_total_element; ++i) {
//     auto el = unlinearized(i);
//     vertex_loads(i) = vertex_load_func(el);

//     UInt nb_adj = dxadj(i + 1) - dxadj(i);
//     for (UInt j = 0; j < nb_adj; ++j, ++adj) {
//       auto el_adj_id = dadjncy(dxadj(i) + j);
//       auto el_adj = unlinearized(el_adj_id);

//       Int load = edge_load_func(el, el_adj);
//       edge_loads(adj) = load;
//     }
//   }

//   AKANTU_DEBUG_OUT();
// }

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

    auto & partition = partitions.alloc(nb_element, 1, type, _not_ghost);
    auto & ghost_part_csr = ghost_partitions_csr(type, _not_ghost);
    ghost_part_csr.resizeRows(nb_element);

    auto & ghost_partition_offset =
        ghost_partitions_offset.alloc(nb_element + 1, 1, type, _ghost);
    auto & ghost_partition = ghost_partitions.alloc(0, 1, type, _ghost);

    const auto & connectivity = mesh.getConnectivity(type, _not_ghost);
    auto conn_it = connectivity.begin(connectivity.getNbComponent());

    for (UInt el = 0; el < nb_element; ++el, ++linearized_el) {
      UInt part = linearized_partitions[linearized_el];

      partition(el) = part;
      std::list<UInt> list_adj_part;
      for (UInt n = 0; n < nb_nodes_per_element; ++n) {
        auto conn = Vector<UInt>(*(conn_it + el));
        UInt node = conn(n);
        for (const auto & adj_element : node_to_elem.getRow(node)) {
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

        ghost_partition.push_back(adj_part);
        ghost_partition_offset(el)++;
      }
    }

    ghost_part_csr.countToCSR();

    /// convert the ghost_partitions_offset array in an offset array
    auto & ghost_partitions_offset_ptr = ghost_partitions_offset(type, _ghost);
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

      auto & partition = partitions.alloc(nb_element, 1, type, _not_ghost);
      AKANTU_DEBUG_INFO("Allocating partitions for " << type);
      auto & ghost_part_csr = ghost_partitions_csr(type, _not_ghost);
      ghost_part_csr.resizeRows(nb_element);

      auto & ghost_partition_offset =
          ghost_partitions_offset.alloc(nb_element + 1, 1, type, _ghost);
      auto & ghost_partition = ghost_partitions.alloc(0, 1, type, _ghost);
      AKANTU_DEBUG_INFO("Allocating ghost_partitions for " << type);
      const Array<std::vector<Element>> & elem_to_subelem =
          mesh.getElementToSubelement(type, _not_ghost);

      // Facet loop
      for (UInt i(0); i < mesh.getNbElement(type, _not_ghost); ++i) {
        const auto & adjacent_elems = elem_to_subelem(i);
        if (adjacent_elems.empty()) {
          partition(i) = 0;
          continue;
        }
        Element min_elem{_max_element_type, std::numeric_limits<UInt>::max(),
              *(ghost_type_t{}.end())};
        UInt min_part(std::numeric_limits<UInt>::max());
        std::set<UInt> adjacent_parts;

        for (auto adj_elem : adjacent_elems) {
          if (adj_elem == ElementNull) { // case of boundary elements
            continue;
          }

          auto adjacent_elem_part = partitions(adj_elem);
          if (adjacent_elem_part < min_part) {
            min_part = adjacent_elem_part;
            min_elem = adj_elem;
          }
          adjacent_parts.insert(adjacent_elem_part);
        }
        partition(i) = min_part;

        auto git = ghost_partitions_csr(min_elem.type, _not_ghost)
                       .begin(min_elem.element);
        auto gend = ghost_partitions_csr(min_elem.type, _not_ghost)
                        .end(min_elem.element);
        for (; git != gend; ++git) {
          adjacent_parts.insert(*git);
        }

        adjacent_parts.erase(min_part);
        for (auto & part : adjacent_parts) {
          ghost_part_csr.getRows().push_back(part);
          ghost_part_csr.rowOffset(i)++;
          ghost_partition.push_back(part);
        }

        ghost_partition_offset(i + 1) =
            ghost_partition_offset(i + 1) + adjacent_elems.size();
      }
      ghost_part_csr.countToCSR();
    }
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshPartition::tweakConnectivity() {
  AKANTU_DEBUG_IN();

  MeshAccessor mesh_accessor(const_cast<Mesh &>(mesh));

  for (auto && type :
       mesh.elementTypes(spatial_dimension, _not_ghost, _ek_not_defined)) {
    auto & connectivity = mesh_accessor.getConnectivity(type, _not_ghost);

    auto & saved_conn = saved_connectivity.alloc(
        connectivity.size(), connectivity.getNbComponent(), type, _not_ghost);
    saved_conn.copy(connectivity);

    for (auto && conn :
         make_view(connectivity, connectivity.getNbComponent())) {
      for (auto && node : conn) {
        if (mesh.isPeriodicSlave(node)) {
          node = mesh.getPeriodicMaster(node);
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
void MeshPartition::printself(std::ostream & stream, int indent) const {
  std::string space(indent, AKANTU_INDENT);
  stream << space << "MeshPartition ["
         << "\n";
  stream << space << " + id           : " << id << "\n";
  stream << space << " + nb partitions: " << nb_partitions << "\n";
  stream << space << " + partitions [ "
         << "\n";
  partitions.printself(stream, indent + 2);
  stream << space << " ]"
         << "\n";
  stream << space << "]"
         << "\n";
}

/* -------------------------------------------------------------------------- */

} // namespace akantu
