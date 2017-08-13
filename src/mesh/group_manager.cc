/**
 * @file   group_manager.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Dana Christen <dana.christen@gmail.com>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Wed Nov 13 2013
 * @date last modification: Mon Aug 17 2015
 *
 * @brief  Stores information about ElementGroup and NodeGroup
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
#include "group_manager.hh"
#include "aka_csr.hh"
#include "data_accessor.hh"
#include "element_group.hh"
#include "element_synchronizer.hh"
#include "mesh.hh"
#include "mesh_accessor.hh"
#include "mesh_utils.hh"
#include "node_group.hh"
/* -------------------------------------------------------------------------- */
#include <algorithm>
#include <iterator>
#include <list>
#include <numeric>
#include <queue>
#include <sstream>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
GroupManager::GroupManager(const Mesh & mesh, const ID & id,
                           const MemoryID & mem_id)
    : id(id), memory_id(mem_id), mesh(mesh) {

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
GroupManager::~GroupManager() {
  auto eit = element_groups.begin();
  auto eend = element_groups.end();
  for (; eit != eend; ++eit)
    delete (eit->second);

  auto nit = node_groups.begin();
  auto nend = node_groups.end();
  for (; nit != nend; ++nit)
    delete (nit->second);
}

/* -------------------------------------------------------------------------- */
NodeGroup & GroupManager::createNodeGroup(const std::string & group_name,
                                          bool replace_group) {
  AKANTU_DEBUG_IN();

  auto it = node_groups.find(group_name);

  if (it != node_groups.end()) {
    if (replace_group) {
      it->second->empty();
      AKANTU_DEBUG_OUT();
      return *(it->second);
    } else
      AKANTU_EXCEPTION(
          "Trying to create a node group that already exists:" << group_name);
  }

  std::stringstream sstr;
  sstr << this->id << ":" << group_name << "_node_group";

  NodeGroup * node_group =
      new NodeGroup(group_name, mesh, sstr.str(), memory_id);

  node_groups[group_name] = node_group;

  AKANTU_DEBUG_OUT();

  return *node_group;
}

/* -------------------------------------------------------------------------- */
template <typename T>
NodeGroup &
GroupManager::createFilteredNodeGroup(const std::string & group_name,
                                      const NodeGroup & source_node_group,
                                      T & filter) {
  AKANTU_DEBUG_IN();

  NodeGroup & node_group = this->createNodeGroup(group_name);
  node_group.append(source_node_group);
  if (T::type == FilterFunctor::_node_filter_functor) {
    node_group.applyNodeFilter(filter);
  } else {
    AKANTU_DEBUG_ERROR("ElementFilter cannot be applied to NodeGroup yet."
                       << " Needs to be implemented.");
  }

  AKANTU_DEBUG_OUT();
  return node_group;
}

/* -------------------------------------------------------------------------- */
void GroupManager::destroyNodeGroup(const std::string & group_name) {
  AKANTU_DEBUG_IN();

  NodeGroups::iterator nit = node_groups.find(group_name);
  NodeGroups::iterator nend = node_groups.end();
  if (nit != nend) {
    delete (nit->second);
    node_groups.erase(nit);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
ElementGroup & GroupManager::createElementGroup(const std::string & group_name,
                                                UInt dimension,
                                                bool replace_group) {
  AKANTU_DEBUG_IN();

  NodeGroup & new_node_group =
      createNodeGroup(group_name + "_nodes", replace_group);

  auto it = element_groups.find(group_name);

  if (it != element_groups.end()) {
    if (replace_group) {
      it->second->empty();
      AKANTU_DEBUG_OUT();
      return *(it->second);
    } else
      AKANTU_EXCEPTION("Trying to create a element group that already exists:"
                       << group_name);
  }

  std::stringstream sstr;
  sstr << this->id << ":" << group_name << "_element_group";

  ElementGroup * element_group = new ElementGroup(
      group_name, mesh, new_node_group, dimension, sstr.str(), memory_id);

  std::stringstream sstr_nodes;
  sstr_nodes << group_name << "_nodes";

  node_groups[sstr_nodes.str()] = &new_node_group;
  element_groups[group_name] = element_group;

  AKANTU_DEBUG_OUT();

  return *element_group;
}

/* -------------------------------------------------------------------------- */
void GroupManager::destroyElementGroup(const std::string & group_name,
                                       bool destroy_node_group) {
  AKANTU_DEBUG_IN();

  auto eit = element_groups.find(group_name);
  auto eend = element_groups.end();
  if (eit != eend) {
    if (destroy_node_group)
      destroyNodeGroup(eit->second->getNodeGroup().getName());
    delete (eit->second);
    element_groups.erase(eit);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void GroupManager::destroyAllElementGroups(bool destroy_node_groups) {
  AKANTU_DEBUG_IN();

  auto eit = element_groups.begin();
  auto eend = element_groups.end();
  for (; eit != eend; ++eit) {
    if (destroy_node_groups)
      destroyNodeGroup(eit->second->getNodeGroup().getName());
    delete (eit->second);
  }
  element_groups.clear();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
ElementGroup & GroupManager::createElementGroup(const std::string & group_name,
                                                UInt dimension,
                                                NodeGroup & node_group) {
  AKANTU_DEBUG_IN();

  if (element_groups.find(group_name) != element_groups.end())
    AKANTU_EXCEPTION(
        "Trying to create a element group that already exists:" << group_name);

  ElementGroup * element_group =
      new ElementGroup(group_name, mesh, node_group, dimension,
                       id + ":" + group_name + "_element_group", memory_id);

  element_groups[group_name] = element_group;

  AKANTU_DEBUG_OUT();

  return *element_group;
}

/* -------------------------------------------------------------------------- */
template <typename T>
ElementGroup & GroupManager::createFilteredElementGroup(
    const std::string & group_name, UInt dimension,
    const NodeGroup & node_group, T & filter) {
  AKANTU_DEBUG_IN();

  ElementGroup * element_group = NULL;

  if (T::type == FilterFunctor::_node_filter_functor) {
    NodeGroup & filtered_node_group = this->createFilteredNodeGroup(
        group_name + "_nodes", node_group, filter);
    element_group =
        &(this->createElementGroup(group_name, dimension, filtered_node_group));
  } else if (T::type == FilterFunctor::_element_filter_functor) {
    AKANTU_DEBUG_ERROR(
        "Cannot handle an ElementFilter yet. Needs to be implemented.");
  }

  AKANTU_DEBUG_OUT();

  return *element_group;
}

/* -------------------------------------------------------------------------- */
class ClusterSynchronizer : public DataAccessor<Element> {
  using DistantIDs = std::set<std::pair<UInt, UInt>>;

public:
  ClusterSynchronizer(GroupManager & group_manager, UInt element_dimension,
                      std::string cluster_name_prefix,
                      ElementTypeMapArray<UInt> & element_to_fragment,
                      const ElementSynchronizer & element_synchronizer,
                      UInt nb_cluster)
      : group_manager(group_manager), element_dimension(element_dimension),
        cluster_name_prefix(cluster_name_prefix),
        element_to_fragment(element_to_fragment),
        element_synchronizer(element_synchronizer), nb_cluster(nb_cluster) {}

  UInt synchronize() {
    StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
    UInt rank = comm.whoAmI();
    UInt nb_proc = comm.getNbProc();

    /// find starting index to renumber local clusters
    Array<UInt> nb_cluster_per_proc(nb_proc);
    nb_cluster_per_proc(rank) = nb_cluster;
    comm.allGather(nb_cluster_per_proc);

    starting_index = std::accumulate(nb_cluster_per_proc.begin(),
                                     nb_cluster_per_proc.begin() + rank, 0);

    UInt global_nb_fragment =
        std::accumulate(nb_cluster_per_proc.begin() + rank,
                        nb_cluster_per_proc.end(), starting_index);

    /// create the local to distant cluster pairs with neighbors
    element_synchronizer.synchronizeOnce(*this, _gst_gm_clusters);

    /// count total number of pairs
    Array<int> nb_pairs(nb_proc); // This is potentially a bug for more than
    // 2**31 pairs, but due to a all gatherv after
    // it must be int to match MPI interfaces
    nb_pairs(rank) = distant_ids.size();
    comm.allGather(nb_pairs);

    UInt total_nb_pairs = std::accumulate(nb_pairs.begin(), nb_pairs.end(), 0);

    /// generate pairs global array
    UInt local_pair_index =
        std::accumulate(nb_pairs.storage(), nb_pairs.storage() + rank, 0);

    Array<UInt> total_pairs(total_nb_pairs, 2);

    for (auto & ids : distant_ids) {
      total_pairs(local_pair_index, 0) = ids.first;
      total_pairs(local_pair_index, 1) = ids.second;
      ++local_pair_index;
    }

    /// communicate pairs to all processors
    nb_pairs *= 2;
    comm.allGatherV(total_pairs, nb_pairs);

    /// renumber clusters

    /// generate fragment list
    std::vector<std::set<UInt>> global_clusters;
    UInt total_nb_cluster = 0;

    Array<bool> is_fragment_in_cluster(global_nb_fragment, 1, false);
    std::queue<UInt> fragment_check_list;

    while (total_pairs.size() != 0) {
      /// create a new cluster
      ++total_nb_cluster;
      global_clusters.resize(total_nb_cluster);
      std::set<UInt> & current_cluster = global_clusters[total_nb_cluster - 1];

      UInt first_fragment = total_pairs(0, 0);
      UInt second_fragment = total_pairs(0, 1);
      total_pairs.erase(0);

      fragment_check_list.push(first_fragment);
      fragment_check_list.push(second_fragment);

      while (!fragment_check_list.empty()) {
        UInt current_fragment = fragment_check_list.front();

        UInt * total_pairs_end =
            total_pairs.storage() + total_pairs.size() * 2;

        UInt * fragment_found =
            std::find(total_pairs.storage(), total_pairs_end, current_fragment);

        if (fragment_found != total_pairs_end) {
          UInt position = fragment_found - total_pairs.storage();
          UInt pair = position / 2;
          UInt other_index = (position + 1) % 2;
          fragment_check_list.push(total_pairs(pair, other_index));
          total_pairs.erase(pair);
        } else {
          fragment_check_list.pop();
          current_cluster.insert(current_fragment);
          is_fragment_in_cluster(current_fragment) = true;
        }
      }
    }

    /// add to FragmentToCluster all local fragments
    for (UInt c = 0; c < global_nb_fragment; ++c) {
      if (!is_fragment_in_cluster(c)) {
        ++total_nb_cluster;
        global_clusters.resize(total_nb_cluster);
        std::set<UInt> & current_cluster =
            global_clusters[total_nb_cluster - 1];

        current_cluster.insert(c);
      }
    }

    /// reorganize element groups to match global clusters
    for (UInt c = 0; c < global_clusters.size(); ++c) {

      /// create new element group corresponding to current cluster
      std::stringstream sstr;
      sstr << cluster_name_prefix << "_" << c;
      ElementGroup & cluster =
          group_manager.createElementGroup(sstr.str(), element_dimension, true);

      std::set<UInt>::iterator it = global_clusters[c].begin();
      std::set<UInt>::iterator end = global_clusters[c].end();

      /// append to current element group all fragments that belong to
      /// the same cluster if they exist
      for (; it != end; ++it) {
        Int local_index = *it - starting_index;

        if (local_index < 0 || local_index >= Int(nb_cluster))
          continue;

        std::stringstream tmp_sstr;
        tmp_sstr << "tmp_" << cluster_name_prefix << "_" << local_index;
        auto eg_it = group_manager.element_group_find(tmp_sstr.str());

        AKANTU_DEBUG_ASSERT(eg_it != group_manager.element_group_end(),
                            "Temporary fragment \"" << tmp_sstr.str()
                                                    << "\" not found");

        cluster.append(*(eg_it->second));
        group_manager.destroyElementGroup(tmp_sstr.str(), true);
      }
    }

    return total_nb_cluster;
  }

private:
  /// functions for parallel communications
  inline UInt getNbData(const Array<Element> & elements,
                        const SynchronizationTag & tag) const {
    if (tag == _gst_gm_clusters)
      return elements.size() * sizeof(UInt);

    return 0;
  }

  inline void packData(CommunicationBuffer & buffer,
                       const Array<Element> & elements,
                       const SynchronizationTag & tag) const {
    if (tag != _gst_gm_clusters)
      return;

    Array<Element>::const_iterator<> el_it = elements.begin();
    Array<Element>::const_iterator<> el_end = elements.end();

    for (; el_it != el_end; ++el_it) {

      const Element & el = *el_it;

      /// for each element pack its global cluster index
      buffer << element_to_fragment(el.type, el.ghost_type)(el.element) +
                    starting_index;
    }
  }

  inline void unpackData(CommunicationBuffer & buffer,
                         const Array<Element> & elements,
                         const SynchronizationTag & tag) {
    if (tag != _gst_gm_clusters)
      return;

    Array<Element>::const_iterator<> el_it = elements.begin();
    Array<Element>::const_iterator<> el_end = elements.end();

    for (; el_it != el_end; ++el_it) {
      UInt distant_cluster;

      buffer >> distant_cluster;

      const Element & el = *el_it;
      UInt local_cluster =
          element_to_fragment(el.type, el.ghost_type)(el.element) +
          starting_index;

      distant_ids.insert(std::make_pair(local_cluster, distant_cluster));
    }
  }

private:
  GroupManager & group_manager;
  UInt element_dimension;
  std::string cluster_name_prefix;
  ElementTypeMapArray<UInt> & element_to_fragment;
  const ElementSynchronizer & element_synchronizer;

  UInt nb_cluster;
  DistantIDs distant_ids;

  UInt starting_index;
};

/* -------------------------------------------------------------------------- */
/// \todo this function doesn't work in 1D
UInt GroupManager::createBoundaryGroupFromGeometry() {
  UInt spatial_dimension = mesh.getSpatialDimension();

  return createClusters(spatial_dimension - 1, "boundary");
}

/* -------------------------------------------------------------------------- */
UInt GroupManager::createClusters(
    UInt element_dimension, Mesh & mesh_facets, std::string cluster_name_prefix,
    const GroupManager::ClusteringFilter & filter) {
  return createClusters(element_dimension, cluster_name_prefix, filter, mesh_facets);
}

/* -------------------------------------------------------------------------- */
UInt GroupManager::createClusters(
    UInt element_dimension, std::string cluster_name_prefix,
    const GroupManager::ClusteringFilter & filter) {
  std::unique_ptr<Mesh> mesh_facets;
  if (!mesh_facets && element_dimension > 0) {
    MeshAccessor mesh_accessor(const_cast<Mesh &>(mesh));
    mesh_facets = std::make_unique<Mesh>(mesh.getSpatialDimension(),
                                         mesh_accessor.getNodesSharedPtr(),
                                         "mesh_facets_for_clusters");

    mesh_facets->defineMeshParent(mesh);

    MeshUtils::buildAllFacets(mesh, *mesh_facets, element_dimension,
                              element_dimension - 1);
  }

  return createClusters(element_dimension, cluster_name_prefix, filter, *mesh_facets);
}

/* -------------------------------------------------------------------------- */
//// \todo if needed element list construction can be optimized by
//// templating the filter class
UInt GroupManager::createClusters(UInt element_dimension,
                                  std::string cluster_name_prefix,
                                  const GroupManager::ClusteringFilter & filter,
                                  Mesh & mesh_facets) {
  AKANTU_DEBUG_IN();

  UInt nb_proc = StaticCommunicator::getStaticCommunicator().getNbProc();
  std::string tmp_cluster_name_prefix = cluster_name_prefix;

  ElementTypeMapArray<UInt> * element_to_fragment = nullptr;

  if (nb_proc > 1 && mesh.isDistributed()) {
    element_to_fragment =
        new ElementTypeMapArray<UInt>("element_to_fragment", id, memory_id);

    element_to_fragment->initialize(
        mesh, _nb_component = 1, _spatial_dimension = element_dimension,
        _element_kind = _ek_not_defined, _with_nb_element = true);
    // mesh.initElementTypeMapArray(*element_to_fragment, 1, element_dimension,
    //                              false, _ek_not_defined, true);
    tmp_cluster_name_prefix = "tmp_" + tmp_cluster_name_prefix;
  }

  ElementTypeMapArray<bool> seen_elements("seen_elements", id, memory_id);
  seen_elements.initialize(mesh, _spatial_dimension = element_dimension,
                           _element_kind = _ek_not_defined,
                           _with_nb_element = true);
  // mesh.initElementTypeMapArray(seen_elements, 1, element_dimension, false,
  //                              _ek_not_defined, true);

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {

    GhostType ghost_type = *gt;
    Element el;
    el.ghost_type = ghost_type;

    Mesh::type_iterator type_it =
        mesh.firstType(element_dimension, ghost_type, _ek_not_defined);
    Mesh::type_iterator type_end =
        mesh.lastType(element_dimension, ghost_type, _ek_not_defined);

    for (; type_it != type_end; ++type_it) {
      el.type = *type_it;
      el.kind = Mesh::getKind(*type_it);
      UInt nb_element = mesh.getNbElement(*type_it, ghost_type);
      Array<bool> & seen_elements_array = seen_elements(el.type, ghost_type);

      for (UInt e = 0; e < nb_element; ++e) {
        el.element = e;
        if (!filter(el))
          seen_elements_array(e) = true;
      }
    }
  }

  Array<bool> checked_node(mesh.getNbNodes(), 1, false);

  UInt nb_cluster = 0;

  /// keep looping until all elements are seen
  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {

    GhostType ghost_type = *gt;
    Element uns_el;
    uns_el.ghost_type = ghost_type;

    Mesh::type_iterator type_it =
        mesh.firstType(element_dimension, ghost_type, _ek_not_defined);
    Mesh::type_iterator type_end =
        mesh.lastType(element_dimension, ghost_type, _ek_not_defined);

    for (; type_it != type_end; ++type_it) {
      uns_el.type = *type_it;
      Array<bool> & seen_elements_vec =
          seen_elements(uns_el.type, uns_el.ghost_type);

      for (UInt e = 0; e < seen_elements_vec.size(); ++e) {
        // skip elements that have been already seen
        if (seen_elements_vec(e) == true)
          continue;

        // set current element
        uns_el.element = e;
        seen_elements_vec(e) = true;

        /// create a new cluster
        std::stringstream sstr;
        sstr << tmp_cluster_name_prefix << "_" << nb_cluster;
        ElementGroup & cluster =
            createElementGroup(sstr.str(), element_dimension, true);
        ++nb_cluster;

        // point element are cluster by themself
        if (element_dimension == 0) {
          cluster.add(uns_el);

          UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(uns_el.type);
          Vector<UInt> connect =
              mesh.getConnectivity(uns_el.type, uns_el.ghost_type)
                  .begin(nb_nodes_per_element)[uns_el.element];
          for (UInt n = 0; n < nb_nodes_per_element; ++n) {
            /// add element's nodes to the cluster
            UInt node = connect[n];
            if (!checked_node(node)) {
              cluster.addNode(node);
              checked_node(node) = true;
            }
          }

          continue;
        }

        std::queue<Element> element_to_add;
        element_to_add.push(uns_el);

        /// keep looping until current cluster is complete (no more
        /// connected elements)
        while (!element_to_add.empty()) {

          /// take first element and erase it in the queue
          Element el = element_to_add.front();
          element_to_add.pop();

          /// if parallel, store cluster index per element
          if (nb_proc > 1 && mesh.isDistributed())
            (*element_to_fragment)(el.type, el.ghost_type)(el.element) =
                nb_cluster - 1;

          /// add current element to the cluster
          cluster.add(el);

          const Array<Element> & element_to_facet =
              mesh_facets.getSubelementToElement(el.type, el.ghost_type);

          UInt nb_facet_per_element = element_to_facet.getNbComponent();

          for (UInt f = 0; f < nb_facet_per_element; ++f) {
            const Element & facet = element_to_facet(el.element, f);

            if (facet == ElementNull)
              continue;

            const std::vector<Element> & connected_elements =
                mesh_facets.getElementToSubelement(
                    facet.type, facet.ghost_type)(facet.element);

            for (UInt elem = 0; elem < connected_elements.size(); ++elem) {
              const Element & check_el = connected_elements[elem];

              // check if this element has to be skipped
              if (check_el == ElementNull || check_el == el)
                continue;

              Array<bool> & seen_elements_vec_current =
                  seen_elements(check_el.type, check_el.ghost_type);

              if (seen_elements_vec_current(check_el.element) == false) {
                seen_elements_vec_current(check_el.element) = true;
                element_to_add.push(check_el);
              }
            }
          }

          UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(el.type);
          Vector<UInt> connect = mesh.getConnectivity(el.type, el.ghost_type)
                                     .begin(nb_nodes_per_element)[el.element];
          for (UInt n = 0; n < nb_nodes_per_element; ++n) {
            /// add element's nodes to the cluster
            UInt node = connect[n];
            if (!checked_node(node)) {
              cluster.addNode(node, false);
              checked_node(node) = true;
            }
          }
        }
      }
    }
  }

  if (nb_proc > 1 && mesh.isDistributed()) {
    ClusterSynchronizer cluster_synchronizer(
        *this, element_dimension, cluster_name_prefix, *element_to_fragment,
        this->mesh.getElementSynchronizer(), nb_cluster);
    nb_cluster = cluster_synchronizer.synchronize();
    delete element_to_fragment;
  }

  if (mesh.isDistributed())
    this->synchronizeGroupNames();

  AKANTU_DEBUG_OUT();
  return nb_cluster;
}

/* -------------------------------------------------------------------------- */
template <typename T>
void GroupManager::createGroupsFromMeshData(const std::string & dataset_name) {
  std::set<std::string> group_names;
  const ElementTypeMapArray<T> & datas = mesh.getData<T>(dataset_name);
  typedef typename ElementTypeMapArray<T>::type_iterator type_iterator;

  std::map<std::string, UInt> group_dim;

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {
    type_iterator type_it = datas.firstType(_all_dimensions, *gt);
    type_iterator type_end = datas.lastType(_all_dimensions, *gt);
    for (; type_it != type_end; ++type_it) {
      const Array<T> & dataset = datas(*type_it, *gt);
      UInt nb_element = mesh.getNbElement(*type_it, *gt);
      AKANTU_DEBUG_ASSERT(dataset.size() == nb_element,
                          "Not the same number of elements ("
                              << *type_it << ":" << *gt
                              << ") in the map from MeshData ("
                              << dataset.size() << ") " << dataset_name
                              << " and in the mesh (" << nb_element << ")!");
      for (UInt e(0); e < nb_element; ++e) {
        std::stringstream sstr;
        sstr << dataset(e);
        std::string gname = sstr.str();
        group_names.insert(gname);

        std::map<std::string, UInt>::iterator it = group_dim.find(gname);
        if (it == group_dim.end()) {
          group_dim[gname] = mesh.getSpatialDimension(*type_it);
        } else {
          it->second = std::max(it->second, mesh.getSpatialDimension(*type_it));
        }
      }
    }
  }

  std::set<std::string>::iterator git = group_names.begin();
  std::set<std::string>::iterator gend = group_names.end();
  for (; git != gend; ++git)
    createElementGroup(*git, group_dim[*git]);

  if (mesh.isDistributed())
    this->synchronizeGroupNames();

  Element el;
  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {
    el.ghost_type = *gt;

    type_iterator type_it = datas.firstType(_all_dimensions, *gt);
    type_iterator type_end = datas.lastType(_all_dimensions, *gt);
    for (; type_it != type_end; ++type_it) {
      el.type = *type_it;

      const Array<T> & dataset = datas(*type_it, *gt);
      UInt nb_element = mesh.getNbElement(*type_it, *gt);
      AKANTU_DEBUG_ASSERT(dataset.size() == nb_element,
                          "Not the same number of elements in the map from "
                          "MeshData and in the mesh!");

      UInt nb_nodes_per_element = mesh.getNbNodesPerElement(el.type);

      Array<UInt>::const_iterator<Vector<UInt>> cit =
          mesh.getConnectivity(*type_it, *gt).begin(nb_nodes_per_element);

      for (UInt e(0); e < nb_element; ++e, ++cit) {
        el.element = e;
        std::stringstream sstr;
        sstr << dataset(e);
        ElementGroup & group = getElementGroup(sstr.str());
        group.add(el, false, false);

        const Vector<UInt> & connect = *cit;
        for (UInt n = 0; n < nb_nodes_per_element; ++n) {
          UInt node = connect[n];
          group.addNode(node, false);
        }
      }
    }
  }

  git = group_names.begin();
  for (; git != gend; ++git) {
    getElementGroup(*git).optimize();
  }
}

template void GroupManager::createGroupsFromMeshData<std::string>(
    const std::string & dataset_name);
template void
GroupManager::createGroupsFromMeshData<UInt>(const std::string & dataset_name);

/* -------------------------------------------------------------------------- */
void GroupManager::createElementGroupFromNodeGroup(
    const std::string & name, const std::string & node_group_name,
    UInt dimension) {
  NodeGroup & node_group = getNodeGroup(node_group_name);
  ElementGroup & group = createElementGroup(name, dimension, node_group);

  group.fillFromNodeGroup();
}

/* -------------------------------------------------------------------------- */
void GroupManager::printself(std::ostream & stream, int indent) const {
  std::string space;
  for (Int i = 0; i < indent; i++, space += AKANTU_INDENT)
    ;

  stream << space << "GroupManager [" << std::endl;

  std::set<std::string> node_group_seen;
  for (const_element_group_iterator it(element_group_begin());
       it != element_group_end(); ++it) {
    it->second->printself(stream, indent + 1);
    node_group_seen.insert(it->second->getNodeGroup().getName());
  }

  for (const_node_group_iterator it(node_group_begin()); it != node_group_end();
       ++it) {
    if (node_group_seen.find(it->second->getName()) == node_group_seen.end())
      it->second->printself(stream, indent + 1);
  }

  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
UInt GroupManager::getNbElementGroups(UInt dimension) const {
  if (dimension == _all_dimensions)
    return element_groups.size();

  ElementGroups::const_iterator it = element_groups.begin();
  ElementGroups::const_iterator end = element_groups.end();
  UInt count = 0;
  for (; it != end; ++it)
    count += (it->second->getDimension() == dimension);
  return count;
}

/* -------------------------------------------------------------------------- */
void GroupManager::checkAndAddGroups(CommunicationBuffer & buffer) {
  AKANTU_DEBUG_IN();

  UInt nb_node_group;
  buffer >> nb_node_group;
  AKANTU_DEBUG_INFO("Received " << nb_node_group << " node group names");

  for (UInt ng = 0; ng < nb_node_group; ++ng) {
    std::string node_group_name;
    buffer >> node_group_name;

    if (node_groups.find(node_group_name) == node_groups.end()) {
      this->createNodeGroup(node_group_name);
    }

    AKANTU_DEBUG_INFO("Received node goup name: " << node_group_name);
  }

  UInt nb_element_group;
  buffer >> nb_element_group;
  AKANTU_DEBUG_INFO("Received " << nb_element_group << " element group names");
  for (UInt eg = 0; eg < nb_element_group; ++eg) {
    std::string element_group_name;
    buffer >> element_group_name;
    std::string node_group_name;
    buffer >> node_group_name;
    UInt dim;
    buffer >> dim;

    AKANTU_DEBUG_INFO("Received element group name: "
                      << element_group_name << " corresponding to a "
                      << Int(dim) << "D group with node group "
                      << node_group_name);

    NodeGroup & node_group = *node_groups[node_group_name];

    if (element_groups.find(element_group_name) == element_groups.end()) {
      this->createElementGroup(element_group_name, dim, node_group);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void GroupManager::fillBufferWithGroupNames(
    DynamicCommunicationBuffer & comm_buffer) const {
  AKANTU_DEBUG_IN();

  // packing node group names;
  UInt nb_groups = this->node_groups.size();
  comm_buffer << nb_groups;
  AKANTU_DEBUG_INFO("Sending " << nb_groups << " node group names");

  NodeGroups::const_iterator nnames_it = node_groups.begin();
  NodeGroups::const_iterator nnames_end = node_groups.end();
  for (; nnames_it != nnames_end; ++nnames_it) {
    std::string node_group_name = nnames_it->first;
    comm_buffer << node_group_name;
    AKANTU_DEBUG_INFO("Sending node goupe name: " << node_group_name);
  }

  // packing element group names with there associated node group name
  nb_groups = this->element_groups.size();
  comm_buffer << nb_groups;
  AKANTU_DEBUG_INFO("Sending " << nb_groups << " element group names");
  ElementGroups::const_iterator gnames_it = this->element_groups.begin();
  ElementGroups::const_iterator gnames_end = this->element_groups.end();
  for (; gnames_it != gnames_end; ++gnames_it) {
    ElementGroup & element_group = *(gnames_it->second);
    std::string element_group_name = gnames_it->first;
    std::string node_group_name = element_group.getNodeGroup().getName();
    UInt dim = element_group.getDimension();

    comm_buffer << element_group_name;
    comm_buffer << node_group_name;
    comm_buffer << dim;

    AKANTU_DEBUG_INFO("Sending element group name: "
                      << element_group_name << " corresponding to a "
                      << Int(dim) << "D group with the node group "
                      << node_group_name);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void GroupManager::synchronizeGroupNames() {
  AKANTU_DEBUG_IN();

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  Int nb_proc = comm.getNbProc();
  Int my_rank = comm.whoAmI();

  if (nb_proc == 1)
    return;

  if (my_rank == 0) {
    for (Int p = 1; p < nb_proc; ++p) {
      CommunicationStatus status;
      comm.probe<char>(p, p, status);
      AKANTU_DEBUG_INFO("Got " << printMemorySize<char>(status.size())
                               << " from proc " << p);

      CommunicationBuffer recv_buffer(status.size());
      comm.receive(recv_buffer, p, p);

      this->checkAndAddGroups(recv_buffer);
    }

    DynamicCommunicationBuffer comm_buffer;
    this->fillBufferWithGroupNames(comm_buffer);

    UInt buffer_size = comm_buffer.size();

    comm.broadcast(buffer_size, 0);

    AKANTU_DEBUG_INFO("Initiating broadcast with "
                      << printMemorySize<char>(comm_buffer.size()));
    comm.broadcast(comm_buffer, 0);

  } else {
    DynamicCommunicationBuffer comm_buffer;
    this->fillBufferWithGroupNames(comm_buffer);

    AKANTU_DEBUG_INFO("Sending " << printMemorySize<char>(comm_buffer.size())
                                 << " to proc " << 0);
    comm.send(comm_buffer, 0, my_rank);

    UInt buffer_size = 0;
    comm.broadcast(buffer_size, 0);

    AKANTU_DEBUG_INFO("Receiving broadcast with "
                      << printMemorySize<char>(comm_buffer.size()));
    CommunicationBuffer recv_buffer(buffer_size);
    comm.broadcast(recv_buffer, 0);

    this->checkAndAddGroups(recv_buffer);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
const ElementGroup &
GroupManager::getElementGroup(const std::string & name) const {
  const_element_group_iterator it = element_group_find(name);
  if (it == element_group_end()) {
    AKANTU_EXCEPTION("There are no element groups named "
                     << name << " associated to the group manager: " << id);
  }

  return *(it->second);
}

/* -------------------------------------------------------------------------- */
ElementGroup & GroupManager::getElementGroup(const std::string & name) {
  element_group_iterator it = element_group_find(name);
  if (it == element_group_end()) {
    AKANTU_EXCEPTION("There are no element groups named "
                     << name << " associated to the group manager: " << id);
  }

  return *(it->second);
}

/* -------------------------------------------------------------------------- */
const NodeGroup & GroupManager::getNodeGroup(const std::string & name) const {
  const_node_group_iterator it = node_group_find(name);
  if (it == node_group_end()) {
    AKANTU_EXCEPTION("There are no node groups named "
                     << name << " associated to the group manager: " << id);
  }

  return *(it->second);
}

/* -------------------------------------------------------------------------- */
NodeGroup & GroupManager::getNodeGroup(const std::string & name) {
  node_group_iterator it = node_group_find(name);
  if (it == node_group_end()) {
    AKANTU_EXCEPTION("There are no node groups named "
                     << name << " associated to the group manager: " << id);
  }

  return *(it->second);
}

} // namespace akantu
