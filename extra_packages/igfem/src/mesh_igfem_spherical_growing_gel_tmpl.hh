/**
 * @file   mesh_igfem_spherical_growing_gel.cc
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @date   Mon Sep 21 12:42:11 2015
 *
 * @brief  Implementation of the functions of MeshIgfemSphericalGrowingGel
 *
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "mesh_utils.hh"
#include <algorithm>

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt dim>
MeshIgfemSphericalGrowingGel<dim>::MeshIgfemSphericalGrowingGel(Mesh & mesh)
    : mesh(mesh), nb_nodes_fem(mesh.getNbNodes()), nb_enriched_nodes(0),
      synchronizer(NULL) {}

/* -------------------------------------------------------------------------- */
template <UInt dim> void MeshIgfemSphericalGrowingGel<dim>::init() {
  nb_nodes_fem = mesh.getNbNodes();

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {
    GhostType ghost_type = *gt;

    Mesh::type_iterator tit = mesh.firstType(dim, ghost_type);
    Mesh::type_iterator tend = mesh.lastType(dim, ghost_type);
    for (; tit != tend;
         ++tit) { // loop to add corresponding IGFEM element types
      if (*tit == _triangle_3) {
        mesh.addConnectivityType(_igfem_triangle_4, ghost_type);
        mesh.addConnectivityType(_igfem_triangle_5, ghost_type);
      } else
        AKANTU_ERROR("Not ready for mesh type " << *tit);
    }

    tit = mesh.firstType(dim, ghost_type, _ek_not_defined);
    tend = mesh.lastType(dim, ghost_type, _ek_not_defined);
    for (; tit != tend; ++tit) {
      AKANTU_BOOST_LIST_SWITCH(INSTANTIATOR, ELEMENT_LIST, *tit);
    }
  }
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
void MeshIgfemSphericalGrowingGel<dim>::computeMeshQueryListIntersectionPoint(
    const std::list<SK::Sphere_3> & query_list) {
  /// store number of currently enriched nodes
  this->nb_enriched_nodes = mesh.getNbNodes() - nb_nodes_fem;
  UInt nb_old_nodes = mesh.getNbNodes();

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {
    GhostType ghost_type = *gt;

    Mesh::type_iterator iit = mesh.firstType(dim, ghost_type, _ek_not_defined);
    Mesh::type_iterator iend = mesh.lastType(dim, ghost_type, _ek_not_defined);
    for (; iit != iend; ++iit) {
      MeshAbstractIntersector<SK::Sphere_3> & intersector =
          *intersectors(*iit, ghost_type);
      intersector.constructData(ghost_type);
      intersector.computeMeshQueryListIntersectionPoint(query_list,
                                                        nb_old_nodes);
      const Array<Real> intersection_points_current_type =
          *(intersector.getIntersectionPoints());
      const Array<UInt> & new_node_per_elem = intersector.getNewNodePerElem();

      /// Send the new node event
      UInt new_points = intersection_points_current_type.getSize();
      StaticCommunicator::getStaticCommunicator().allReduce(&new_points, 1,
                                                            _so_sum);

      if (new_points > 0) {
        Array<Real> & nodes = this->mesh.getNodes();
        UInt nb_node = nodes.getSize();

        Array<Real>::const_vector_iterator intersec_p_it =
            intersection_points_current_type.begin(dim);

        NewIGFEMNodesEvent new_nodes_event;
        for (UInt in = 0; in < intersection_points_current_type.getSize();
             ++in, ++intersec_p_it) {
          nodes.push_back(*intersec_p_it);
          new_nodes_event.getList().push_back(nb_node + in);
        }
        new_nodes_event.setNewNodePerElem(new_node_per_elem);
        new_nodes_event.setType(*iit);
        new_nodes_event.setGhostType(ghost_type);
        this->mesh.sendEvent(new_nodes_event);
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
void MeshIgfemSphericalGrowingGel<dim>::removeAdditionalNodes() {
  AKANTU_DEBUG_IN();

  UInt total_nodes = this->mesh.getNbNodes();
  UInt nb_new_enriched_nodes =
      total_nodes - this->nb_enriched_nodes - this->nb_nodes_fem;
  UInt old_total_nodes = this->nb_nodes_fem + this->nb_enriched_nodes;

  UInt total_new_nodes = nb_new_enriched_nodes;
  StaticCommunicator::getStaticCommunicator().allReduce(&total_new_nodes, 1,
                                                        _so_sum);
  if (total_new_nodes == 0)
    return;

  RemovedNodesEvent remove_nodes(this->mesh);
  Array<UInt> & nodes_removed = remove_nodes.getList();
  Array<UInt> & new_numbering = remove_nodes.getNewNumbering();

  for (UInt nnod = 0; nnod < this->nb_nodes_fem; ++nnod) {
    new_numbering(nnod) = nnod;
  }

  for (UInt nnod = nb_nodes_fem; nnod < old_total_nodes; ++nnod) {
    new_numbering(nnod) = UInt(-1);
    nodes_removed.push_back(nnod);
  }

  for (UInt nnod = 0; nnod < nb_new_enriched_nodes; ++nnod) {
    new_numbering(nnod + old_total_nodes) = this->nb_nodes_fem + nnod;
  }

  if (nb_new_enriched_nodes > 0) {
    for (UInt gt = _not_ghost; gt <= _ghost; ++gt) {
      GhostType ghost_type = (GhostType)gt;

      Mesh::type_iterator it =
          mesh.firstType(_all_dimensions, ghost_type, _ek_not_defined);
      Mesh::type_iterator end =
          mesh.lastType(_all_dimensions, ghost_type, _ek_not_defined);
      for (; it != end; ++it) {
        ElementType type = *it;
        Array<UInt> & connectivity_array =
            mesh.getConnectivity(type, ghost_type);
        UInt nb_nodes_per_element = connectivity_array.getNbComponent();

        Array<UInt>::vector_iterator conn_it =
            connectivity_array.begin(nb_nodes_per_element);
        Array<UInt>::vector_iterator conn_end =
            connectivity_array.end(nb_nodes_per_element);

        for (; conn_it != conn_end; ++conn_it) {
          for (UInt n = 0; n < nb_nodes_per_element; ++n) {
            (*conn_it)(n) = new_numbering((*conn_it)(n));
          }
        }
      }
    }
  }

  this->mesh.sendEvent(remove_nodes);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt dim> void MeshIgfemSphericalGrowingGel<dim>::buildIGFEMMesh() {
  AKANTU_DEBUG_IN();

  NewIGFEMElementsEvent new_elements_event;
  UInt total_new_elements = 0;
  Array<Element> & new_elements_list = new_elements_event.getList();
  Array<Element> & old_elements_list = new_elements_event.getOldElementsList();
  RemovedElementsEvent removed_elements_event(this->mesh);
  ChangedElementsEvent changed_elements_event(this->mesh);

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {
    GhostType ghost_type = *gt;

    mesh.initElementTypeMapArray(removed_elements_event.getNewNumbering(), 1,
                                 dim, ghost_type, false, _ek_not_defined);
    mesh.initElementTypeMapArray(changed_elements_event.getNewNumbering(), 1,
                                 dim, ghost_type, false, _ek_not_defined);

    /// loop over all types in the mesh for a given ghost type
    Mesh::type_iterator iit = mesh.firstType(dim, ghost_type, _ek_not_defined);
    Mesh::type_iterator iend = mesh.lastType(dim, ghost_type, _ek_not_defined);
    for (; iit != iend; ++iit) {
      ElementType type = *iit;

      MeshAbstractIntersector<SK::Sphere_3> & intersector =
          *intersectors(type, ghost_type);
      const Array<UInt> & new_node_per_elem = intersector.getNewNodePerElem();

      UInt n_new_el = 0;
      Array<UInt> & connectivity = this->mesh.getConnectivity(type, ghost_type);

      /// get the connectivities of all types that that may transform
      Array<UInt> & connec_igfem_tri_4 =
          this->mesh.getConnectivity(_igfem_triangle_4, ghost_type);
      Array<UInt> & connec_igfem_tri_5 =
          this->mesh.getConnectivity(_igfem_triangle_5, ghost_type);
      Array<UInt> & connec_tri_3 =
          this->mesh.getConnectivity(_triangle_3, ghost_type);

      /// create elements to store the newly generated elements
      Element el_tri_3(_triangle_3, 0, ghost_type, _ek_regular);
      Element el_igfem_tri_4(_igfem_triangle_4, 0, ghost_type, _ek_igfem);
      Element el_igfem_tri5(_igfem_triangle_5, 0, ghost_type, _ek_igfem);

      Array<UInt> & new_numbering =
          removed_elements_event.getNewNumbering(type, ghost_type);
      new_numbering.resize(connectivity.getSize());
      /// container for element to be removed
      Element old_element(type, 0, ghost_type, Mesh::getKind(type));

      for (UInt nel = 0; nel != new_node_per_elem.getSize(); ++nel) {
        /// a former IGFEM triangle is transformed into a regular triangle
        if ((type != _triangle_3) && (new_node_per_elem(nel, 0) == 0)) {

          Vector<UInt> connec_new_elem(3);
          connec_new_elem(0) = connectivity(nel, 0);
          connec_new_elem(1) = connectivity(nel, 1);
          connec_new_elem(2) = connectivity(nel, 2);
          /// add the new element in the mesh
          UInt nb_triangles_3 = connec_tri_3.getSize();
          connec_tri_3.push_back(connec_new_elem);
          el_tri_3.element = nb_triangles_3;
          new_elements_list.push_back(el_tri_3);
          /// the number of the old element in the mesh
          old_element.element = nel;
          old_elements_list.push_back(old_element);
          new_numbering(nel) = UInt(-1);
          ++n_new_el;
        }

        /// element is enriched
        else if (new_node_per_elem(nel, 0) != 0) {
          switch (new_node_per_elem(nel, 0)) {
          /// new element is of type igfem_triangle_4
          case 1: {
            Vector<UInt> connec_new_elem(4);
            switch (new_node_per_elem(nel, 2)) {
            case 0:
              connec_new_elem(0) = connectivity(nel, 2);
              connec_new_elem(1) = connectivity(nel, 0);
              connec_new_elem(2) = connectivity(nel, 1);
              break;
            case 1:
              connec_new_elem(0) = connectivity(nel, 0);
              connec_new_elem(1) = connectivity(nel, 1);
              connec_new_elem(2) = connectivity(nel, 2);
              break;
            case 2:
              connec_new_elem(0) = connectivity(nel, 1);
              connec_new_elem(1) = connectivity(nel, 2);
              connec_new_elem(2) = connectivity(nel, 0);
              break;
            default:
              AKANTU_ERROR("A triangle has only 3 segments not "
                           << new_node_per_elem(nel, 2));
              break;
            }
            connec_new_elem(3) = new_node_per_elem(nel, 1);
            UInt nb_igfem_triangles_4 = connec_igfem_tri_4.getSize();
            connec_igfem_tri_4.push_back(connec_new_elem);
            el_igfem_tri_4.element = nb_igfem_triangles_4;
            new_elements_list.push_back(el_igfem_tri_4);
            break;
          }
          /// new element is of type igfem_triangle_5
          case 2: {
            Vector<UInt> connec_new_elem(5);
            if ((new_node_per_elem(nel, 2) == 0) &&
                (new_node_per_elem(nel, 4) == 1)) {
              connec_new_elem(0) = connectivity(nel, 1);
              connec_new_elem(1) = connectivity(nel, 2);
              connec_new_elem(2) = connectivity(nel, 0);
              connec_new_elem(3) = new_node_per_elem(nel, 3);
              connec_new_elem(4) = new_node_per_elem(nel, 1);
            } else if ((new_node_per_elem(nel, 2) == 0) &&
                       (new_node_per_elem(nel, 4) == 2)) {
              connec_new_elem(0) = connectivity(nel, 0);
              connec_new_elem(1) = connectivity(nel, 1);
              connec_new_elem(2) = connectivity(nel, 2);
              connec_new_elem(3) = new_node_per_elem(nel, 1);
              connec_new_elem(4) = new_node_per_elem(nel, 3);
            } else if ((new_node_per_elem(nel, 2) == 1) &&
                       (new_node_per_elem(nel, 4) == 2)) {
              connec_new_elem(0) = connectivity(nel, 2);
              connec_new_elem(1) = connectivity(nel, 0);
              connec_new_elem(2) = connectivity(nel, 1);
              connec_new_elem(3) = new_node_per_elem(nel, 3);
              connec_new_elem(4) = new_node_per_elem(nel, 1);
            } else if ((new_node_per_elem(nel, 2) == 1) &&
                       (new_node_per_elem(nel, 4) == 0)) {
              connec_new_elem(0) = connectivity(nel, 1);
              connec_new_elem(1) = connectivity(nel, 2);
              connec_new_elem(2) = connectivity(nel, 0);
              connec_new_elem(3) = new_node_per_elem(nel, 1);
              connec_new_elem(4) = new_node_per_elem(nel, 3);
            } else if ((new_node_per_elem(nel, 2) == 2) &&
                       (new_node_per_elem(nel, 4) == 0)) {
              connec_new_elem(0) = connectivity(nel, 0);
              connec_new_elem(1) = connectivity(nel, 1);
              connec_new_elem(2) = connectivity(nel, 2);
              connec_new_elem(3) = new_node_per_elem(nel, 3);
              connec_new_elem(4) = new_node_per_elem(nel, 1);
            } else if ((new_node_per_elem(nel, 2) == 2) &&
                       (new_node_per_elem(nel, 4) == 1)) {
              connec_new_elem(0) = connectivity(nel, 2);
              connec_new_elem(1) = connectivity(nel, 0);
              connec_new_elem(2) = connectivity(nel, 1);
              connec_new_elem(3) = new_node_per_elem(nel, 1);
              connec_new_elem(4) = new_node_per_elem(nel, 3);
            } else
              AKANTU_ERROR("A triangle has only 3 segments (0 to 2) not "
                           << new_node_per_elem(nel, 2) << "and"
                           << new_node_per_elem(nel, 4));

            UInt nb_igfem_triangles_5 = connec_igfem_tri_5.getSize();
            connec_igfem_tri_5.push_back(connec_new_elem);
            el_igfem_tri5.element = nb_igfem_triangles_5;
            new_elements_list.push_back(el_igfem_tri5);
            break;
          }
          default:
            AKANTU_ERROR("IGFEM cannot add " << new_node_per_elem(nel, 0)
                                             << " nodes to triangles");
            break;
          }
          old_element.element = nel;
          old_elements_list.push_back(old_element);
          new_numbering(nel) = UInt(-1);
          ++n_new_el;
        }
      }
      total_new_elements += n_new_el;
    }
  }

  /// update new numbering
  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {
    GhostType ghost_type = *gt;
    /// loop over all types in the mesh for a given ghost type
    Mesh::type_iterator iit = mesh.firstType(dim, ghost_type, _ek_not_defined);
    Mesh::type_iterator iend = mesh.lastType(dim, ghost_type, _ek_not_defined);
    for (; iit != iend; ++iit) {
      ElementType type = *iit;
      Array<UInt> & new_numbering =
          removed_elements_event.getNewNumbering(type, ghost_type);
      UInt el_index = 0;
      UInt nb_element = this->mesh.getNbElement(type, ghost_type);
      new_numbering.resize(nb_element);
      for (UInt e = 0; e < nb_element; ++e) {
        if (new_numbering(e) != UInt(-1)) {
          new_numbering(e) = el_index;
          ++el_index;
        }
      }
      changed_elements_event.getNewNumbering(type, ghost_type)
          .copy(new_numbering);
    }
  }

  StaticCommunicator::getStaticCommunicator().allReduce(&total_new_elements, 1,
                                                        _so_sum);

  if (total_new_elements > 0) {
    changed_elements_event.getListOld().copy(
        new_elements_event.getOldElementsList());
    changed_elements_event.getListNew().copy(new_elements_event.getList());
    this->mesh.sendEvent(changed_elements_event);

    this->mesh.sendEvent(new_elements_event);

    Array<Element> & removed_list = removed_elements_event.getList();
    removed_list.copy(new_elements_event.getOldElementsList());
    this->mesh.sendEvent(removed_elements_event);
  }

  removeAdditionalNodes();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
void MeshIgfemSphericalGrowingGel<dim>::buildSegmentConnectivityToNodeType() {
  Mesh mesh_facets(mesh.initMeshFacets());
  MeshUtils::buildSegmentToNodeType(mesh, mesh_facets, synchronizer);

  // only the ghost elements are considered
  for (UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType ghost_type = (GhostType)g;

    Mesh::type_iterator it = mesh_facets.firstType(1, ghost_type);
    Mesh::type_iterator end = mesh_facets.lastType(1, ghost_type);
    for (; it != end; ++it) {
      ElementType type = *it;

      const Array<Int> & segment_to_nodetype =
          mesh_facets.getData<Int>("segment_to_nodetype", type, ghost_type);

      const Array<UInt> & segment_connectivity =
          mesh_facets.getConnectivity(type, ghost_type);

      // looping over all the segments
      Array<UInt>::const_vector_iterator conn_it =
          segment_connectivity.begin(segment_connectivity.getNbComponent());
      Array<UInt>::const_vector_iterator conn_end =
          segment_connectivity.end(segment_connectivity.getNbComponent());
      UInt seg_index = 0;

      UInt n[2];
      for (; conn_it != conn_end; ++conn_it, ++seg_index) {
        Int seg_type = segment_to_nodetype(seg_index);
        n[0] = (*conn_it)(0);
        n[1] = (*conn_it)(1);

        if ((mesh.isMasterNode(n[0]) || mesh.isSlaveNode(n[0])) &&
            (mesh.isMasterNode(n[1]) || mesh.isSlaveNode(n[1]))) {
          std::sort(n, n + 2);
          segment_conn_to_node_type[std::make_pair(n[0], n[1])] = seg_type;
        }
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
void MeshIgfemSphericalGrowingGel<dim>::updateNodeType(
    const Array<UInt> & nodes_list, const Array<UInt> & new_node_per_elem,
    ElementType type, GhostType ghost_type) {
  Array<Int> & nodes_type = mesh.getNodesType();
  UInt old_nodes = nodes_type.getSize();
  UInt new_nodes = nodes_list.getSize();

  // exit this function if the simulation in run in serial
  if (old_nodes == 0 || new_nodes == 0)
    return;

  nodes_type.resize(old_nodes + new_nodes);
  Array<bool> checked_node(new_nodes, 1, false);

  UInt nb_prim_by_el = 0;
  if ((type == _triangle_3) || (type == _igfem_triangle_4) ||
      (type == _igfem_triangle_5)) {
    nb_prim_by_el = 3;
  } else {
    AKANTU_ERROR("Not ready for mesh type " << type);
  }

  // determine the node type for the local, master and slave nodes
  const Array<UInt> & connectivity = mesh.getConnectivity(type, ghost_type);
  unordered_map<std::pair<UInt, UInt>, Int>::type::iterator seg_type_it;
  unordered_map<std::pair<UInt, UInt>, Int>::type::iterator seg_type_end =
      segment_conn_to_node_type.end();

  for (UInt el = 0; el < new_node_per_elem.getSize(); ++el) {
    UInt nb_nodes = new_node_per_elem(el, 0);

    for (UInt n = 0; n < nb_nodes; ++n) {
      UInt node_index = new_node_per_elem(el, 2 * n + 1);
      if (node_index < old_nodes || checked_node(node_index - old_nodes))
        continue;

      // get the elemental connectivity of the segment associated to the node
      UInt segment_index = new_node_per_elem(el, 2 * n + 2);

      UInt extreme_nodes[2];
      extreme_nodes[0] = segment_index;
      extreme_nodes[1] = segment_index + 1;
      if (extreme_nodes[1] == nb_prim_by_el)
        extreme_nodes[1] = 0;

      // get the connectivity of the segment
      extreme_nodes[0] = connectivity(el, extreme_nodes[0]);
      extreme_nodes[1] = connectivity(el, extreme_nodes[1]);

      // if one extreme nodes is pure ghost, then also the new node is pure
      // ghost
      if (mesh.isPureGhostNode(extreme_nodes[0]) ||
          mesh.isPureGhostNode(extreme_nodes[1]))
        nodes_type(node_index) = -3;
      // if on of the two extreme nodes is local, then also the new node is
      // local
      else if (mesh.isLocalNode(extreme_nodes[0]) ||
               mesh.isLocalNode(extreme_nodes[1]))
        nodes_type(node_index) = -1;

      // otherwise use the value stored in the map
      else {
        std::sort(extreme_nodes, extreme_nodes + 2);

        seg_type_it = segment_conn_to_node_type.find(
            std::make_pair(extreme_nodes[0], extreme_nodes[1]));

        AKANTU_DEBUG_ASSERT(seg_type_it != seg_type_end, "Segment not found");

        nodes_type(node_index) = seg_type_it->second;
      }

      checked_node(node_index - old_nodes) = true;
    }
  }

  AKANTU_DEBUG_ASSERT(std::accumulate(checked_node.begin(), checked_node.end(),
                                      0) == Int(checked_node.getSize()),
                      "Not all new nodes were assigned a node type");
}

} // namespace akantu
