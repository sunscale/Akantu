/**
 * @file   grid_synchronizer.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Oct 03 2011
 * @date last modification: Tue Nov 07 2017
 *
 * @brief  implementation of the grid synchronizer
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
#include "grid_synchronizer.hh"
#include "aka_grid_dynamic.hh"
#include "communicator.hh"
#include "fe_engine.hh"
#include "integration_point.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include <iostream>

/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <class E>
void GridSynchronizer::createGridSynchronizer(const SpatialGrid<E> & grid) {
  AKANTU_DEBUG_IN();

  const Communicator & comm = this->mesh.getCommunicator();
  UInt nb_proc = comm.getNbProc();
  UInt my_rank = comm.whoAmI();

  if (nb_proc == 1) {
    return;
  }

  UInt spatial_dimension = this->mesh.getSpatialDimension();

  BBox my_bounding_box(spatial_dimension);

  const auto & lower = grid.getLowerBounds();
  const auto & upper = grid.getUpperBounds();
  const auto & spacing = grid.getSpacing();

  my_bounding_box.getLowerBounds() = lower - spacing;
  my_bounding_box.getUpperBounds() = upper + spacing;

  AKANTU_DEBUG_INFO(
      "Exchange of bounding box to detect the overlapping regions.");

  auto && bboxes = my_bounding_box.allGather(comm);

  std::vector<bool> intersects_proc(nb_proc);
  std::fill(intersects_proc.begin(), intersects_proc.end(), true);

  Matrix<Int> first_cells(spatial_dimension, nb_proc);
  Matrix<Int> last_cells(spatial_dimension, nb_proc);

  std::map<UInt, ElementTypeMapArray<UInt>> element_per_proc;

  // check the overlapping between my box and the one from other processors
  for (UInt p = 0; p < nb_proc; ++p) {
    if (p == my_rank) {
      continue;
    }

    const auto & proc_bounding_box = bboxes[p];
    auto intersection = my_bounding_box.intersection(proc_bounding_box);

    Vector<Int> first_cell_p = first_cells(p);
    Vector<Int> last_cell_p = last_cells(p);

    intersects_proc[p] = intersection;

    if (intersects_proc[p]) {
      for (UInt s = 0; s < spatial_dimension; ++s) {
        first_cell_p(s) = grid.getCellID(intersection.getLowerBounds()(s), s);
        last_cell_p(s) = grid.getCellID(intersection.getUpperBounds()(s), s);
      }
    }

    // create the list of cells in the overlapping
    using CellID = typename SpatialGrid<E>::CellID;

    std::vector<CellID> cell_ids;

    if (intersects_proc[p]) {
      AKANTU_DEBUG_INFO("I intersects with processor " << p);

      CellID cell_id(spatial_dimension);

      // for (UInt i = 0; i < spatial_dimension; ++i) {
      //   if(first_cell_p[i] != 0) --first_cell_p[i];
      //   if(last_cell_p[i] != 0) ++last_cell_p[i];
      // }

      for (Int fd = first_cell_p(0); fd <= last_cell_p(0); ++fd) {
        cell_id.setID(0, fd);
        if (spatial_dimension == 1) {
          cell_ids.push_back(cell_id);
        } else {
          for (Int sd = first_cell_p(1); sd <= last_cell_p(1); ++sd) {
            cell_id.setID(1, sd);
            if (spatial_dimension == 2) {
              cell_ids.push_back(cell_id);
            } else {
              for (Int ld = first_cell_p(2); ld <= last_cell_p(2); ++ld) {
                cell_id.setID(2, ld);
                cell_ids.push_back(cell_id);
              }
            }
          }
        }
      }

      // get the list of elements in the cells of the overlapping
      std::set<Element> to_send;
      for (auto & cur_cell_id : cell_ids) {
        auto & cell = grid.getCell(cur_cell_id);
        for (auto & element : cell) {
          to_send.insert(element);
        }
      }

      AKANTU_DEBUG_INFO("I have prepared " << to_send.size()
                                           << " elements to send to processor "
                                           << p);
      auto & scheme = this->getCommunications().createSendScheme(p);
      std::stringstream sstr;
      sstr << "element_per_proc_" << p;
      element_per_proc.emplace(
          std::piecewise_construct, std::forward_as_tuple(p),
          std::forward_as_tuple(sstr.str(), id));

      ElementTypeMapArray<UInt> & elempproc = element_per_proc[p];

      for (auto elem : to_send) {
        ElementType type = elem.type;
        UInt nb_nodes_per_element = mesh.getNbNodesPerElement(type);

        // /!\ this part must be slow due to the access in the
        // ElementTypeMapArray<UInt>
        if (!elempproc.exists(type, _not_ghost)) {
          elempproc.alloc(0, nb_nodes_per_element, type, _not_ghost);
        }

        Vector<UInt> global_connect(nb_nodes_per_element);
        Vector<UInt> local_connect = mesh.getConnectivity(type).begin(
            nb_nodes_per_element)[elem.element];

        for (UInt i = 0; i < nb_nodes_per_element; ++i) {
          global_connect(i) = mesh.getNodeGlobalId(local_connect(i));
          AKANTU_DEBUG_ASSERT(
              global_connect(i) < mesh.getNbGlobalNodes(),
              "This global node send in the connectivity does not seem correct "
                  << global_connect(i) << " corresponding to "
                  << local_connect(i) << " from element " << elem.element);
        }

        elempproc(type).push_back(global_connect);
        scheme.push_back(elem);
      }
    }
  }

  AKANTU_DEBUG_INFO("I have finished to compute intersection,"
                    << " no it's time to communicate with my neighbors");

  /**
   * Sending loop, sends the connectivity asynchronously to all concerned proc
   */
  std::vector<CommunicationRequest> isend_requests;
  Tensor3<UInt> space(2, _max_element_type, nb_proc);

  for (UInt p = 0; p < nb_proc; ++p) {
    if (p == my_rank) {
      continue;
    }

    if (not intersects_proc[p]) {
      continue;
    }

    Matrix<UInt> info_proc = space(p);
    auto & elempproc = element_per_proc[p];
    UInt count = 0;

    for (auto type : elempproc.elementTypes(_all_dimensions, _not_ghost)) {
      Array<UInt> & conn = elempproc(type, _not_ghost);

      Vector<UInt> info = info_proc((UInt)type);
      info[0] = (UInt)type;
      info[1] = conn.size() * conn.getNbComponent();

      AKANTU_DEBUG_INFO(
          "I have " << conn.size() << " elements of type " << type
                    << " to send to processor " << p << " (communication tag : "
                    << Tag::genTag(my_rank, count, DATA_TAG) << ")");

      isend_requests.push_back(
          comm.asyncSend(info, p, Tag::genTag(my_rank, count, SIZE_TAG)));

      if (info[1] != 0) {
        isend_requests.push_back(comm.asyncSend<UInt>(
            conn, p, Tag::genTag(my_rank, count, DATA_TAG)));
      }

      ++count;
    }

    Vector<UInt> info = info_proc((UInt)_not_defined);
    info[0] = (UInt)_not_defined;
    info[1] = 0;
    isend_requests.push_back(
        comm.asyncSend(info, p, Tag::genTag(my_rank, count, SIZE_TAG)));
  }

  /**
   * Receives the connectivity and store them in the ghosts elements
   */
  MeshAccessor mesh_accessor(mesh);
  auto & global_nodes_ids = mesh_accessor.getNodesGlobalIds();
  auto & nodes_type = mesh_accessor.getNodesFlags();

  std::vector<CommunicationRequest> isend_nodes_requests;
  Vector<UInt> nb_nodes_to_recv(nb_proc);

  UInt nb_total_nodes_to_recv = 0;
  UInt nb_current_nodes = global_nodes_ids.size();

  NewNodesEvent new_nodes;
  NewElementsEvent new_elements;

  std::map<UInt, std::vector<UInt>> ask_nodes_per_proc;

  for (UInt p = 0; p < nb_proc; ++p) {
    nb_nodes_to_recv(p) = 0;
    if (p == my_rank) {
      continue;
    }

    if (!intersects_proc[p]) {
      continue;
    }

    auto & scheme = this->getCommunications().createRecvScheme(p);

    ask_nodes_per_proc.emplace(std::piecewise_construct,
                               std::forward_as_tuple(p),
                               std::forward_as_tuple(0));

    auto & ask_nodes = ask_nodes_per_proc[p];
    UInt count = 0;

    ElementType type = _not_defined;
    do {
      Vector<UInt> info(2);
      comm.receive(info, p, Tag::genTag(p, count, SIZE_TAG));

      type = (ElementType)info[0];

      if (type == _not_defined) {
        break;
      }

      UInt nb_nodes_per_element = mesh.getNbNodesPerElement(type);
      UInt nb_element = info[1] / nb_nodes_per_element;

      Array<UInt> tmp_conn(nb_element, nb_nodes_per_element);
      tmp_conn.zero();
      if (info[1] != 0) {
        comm.receive<UInt>(tmp_conn, p, Tag::genTag(p, count, DATA_TAG));
      }

      AKANTU_DEBUG_INFO("I will receive "
                        << nb_element << " elements of type "
                        << ElementType(info[0]) << " from processor " << p
                        << " (communication tag : "
                        << Tag::genTag(p, count, DATA_TAG) << ")");

      auto & ghost_connectivity = mesh_accessor.getConnectivity(type, _ghost);
      auto & ghost_counter = mesh_accessor.getGhostsCounters(type, _ghost);

      UInt nb_ghost_element = ghost_connectivity.size();
      Element element{type, 0, _ghost};

      Vector<UInt> conn(nb_nodes_per_element);
      for (UInt el = 0; el < nb_element; ++el) {
        UInt nb_node_to_ask_for_elem = 0;

        for (UInt n = 0; n < nb_nodes_per_element; ++n) {
          UInt gn = tmp_conn(el, n);
          UInt ln = global_nodes_ids.find(gn);

          AKANTU_DEBUG_ASSERT(gn < mesh.getNbGlobalNodes(),
                              "This global node seems not correct "
                                  << gn << " from element " << el << " node "
                                  << n);

          if (ln == UInt(-1)) {
            global_nodes_ids.push_back(gn);
            nodes_type.push_back(NodeFlag::_pure_ghost); // pure ghost node
            ln = nb_current_nodes;

            new_nodes.getList().push_back(ln);
            ++nb_current_nodes;
            ask_nodes.push_back(gn);
            ++nb_node_to_ask_for_elem;
          }
          conn[n] = ln;
        }

        // all the nodes are already known locally, the element should
        // already exists
        auto c = UInt(-1);
        if (nb_node_to_ask_for_elem == 0) {
          c = ghost_connectivity.find(conn);
          element.element = c;
        }

        if (c == UInt(-1)) {
          element.element = nb_ghost_element;
          ++nb_ghost_element;
          ghost_connectivity.push_back(conn);
          ghost_counter.push_back(1);
          new_elements.getList().push_back(element);
        } else {
          ++ghost_counter(c);
        }

        scheme.push_back(element);
      }

      count++;
    } while (type != _not_defined);

    AKANTU_DEBUG_INFO("I have "
                      << ask_nodes.size()
                      << " missing nodes for elements coming from processor "
                      << p << " (communication tag : "
                      << Tag::genTag(my_rank, 0, ASK_NODES_TAG) << ")");

    ask_nodes.push_back(UInt(-1));

    isend_nodes_requests.push_back(
        comm.asyncSend(ask_nodes, p, Tag::genTag(my_rank, 0, ASK_NODES_TAG)));
    nb_nodes_to_recv(p) = ask_nodes.size() - 1;
    nb_total_nodes_to_recv += nb_nodes_to_recv(p);
  }

  Communicator::waitAll(isend_requests);
  Communicator::freeCommunicationRequest(isend_requests);

  /**
   * Sends requested nodes to proc
   */
  auto & nodes = const_cast<Array<Real> &>(mesh.getNodes());
  UInt nb_nodes = nodes.size();

  std::vector<CommunicationRequest> isend_coordinates_requests;
  std::map<UInt, Array<Real>> nodes_to_send_per_proc;
  for (UInt p = 0; p < nb_proc; ++p) {
    if (p == my_rank || !intersects_proc[p]) {
      continue;
    }

    Array<UInt> asked_nodes;
    CommunicationStatus status;
    AKANTU_DEBUG_INFO("Waiting list of nodes to send to processor "
                      << p << "(communication tag : "
                      << Tag::genTag(p, 0, ASK_NODES_TAG) << ")");

    comm.probe<UInt>(p, Tag::genTag(p, 0, ASK_NODES_TAG), status);
    UInt nb_nodes_to_send = status.size();
    asked_nodes.resize(nb_nodes_to_send);

    AKANTU_DEBUG_INFO("I have " << nb_nodes_to_send - 1
                                << " nodes to send to processor " << p
                                << " (communication tag : "
                                << Tag::genTag(p, 0, ASK_NODES_TAG) << ")");

    AKANTU_DEBUG_INFO("Getting list of nodes to send to processor "
                      << p << " (communication tag : "
                      << Tag::genTag(p, 0, ASK_NODES_TAG) << ")");

    comm.receive(asked_nodes, p, Tag::genTag(p, 0, ASK_NODES_TAG));

    nb_nodes_to_send--;
    asked_nodes.resize(nb_nodes_to_send);

    nodes_to_send_per_proc.emplace(std::piecewise_construct,
                                   std::forward_as_tuple(p),
                                   std::forward_as_tuple(0, spatial_dimension));

    auto & nodes_to_send = nodes_to_send_per_proc[p];
    auto node_it = nodes.begin(spatial_dimension);

    for (UInt n = 0; n < nb_nodes_to_send; ++n) {
      UInt ln = global_nodes_ids.find(asked_nodes(n));
      AKANTU_DEBUG_ASSERT(ln != UInt(-1), "The node ["
                                              << asked_nodes(n)
                                              << "] requested by proc " << p
                                              << " was not found locally!");
      nodes_to_send.push_back(node_it + ln);
    }

    if (nb_nodes_to_send != 0) {
      AKANTU_DEBUG_INFO("Sending the "
                        << nb_nodes_to_send << " nodes to processor " << p
                        << " (communication tag : "
                        << Tag::genTag(p, 0, SEND_NODES_TAG) << ")");

      isend_coordinates_requests.push_back(comm.asyncSend(
          nodes_to_send, p, Tag::genTag(my_rank, 0, SEND_NODES_TAG)));
    }
#if not defined(AKANTU_NDEBUG)
    else {
      AKANTU_DEBUG_INFO("No nodes to send to processor " << p);
    }
#endif
  }

  Communicator::waitAll(isend_nodes_requests);
  Communicator::freeCommunicationRequest(isend_nodes_requests);

  nodes.resize(nb_total_nodes_to_recv + nb_nodes);
  for (UInt p = 0; p < nb_proc; ++p) {
    if ((p != my_rank) && (nb_nodes_to_recv(p) > 0)) {
      AKANTU_DEBUG_INFO("Receiving the "
                        << nb_nodes_to_recv(p) << " nodes from processor " << p
                        << " (communication tag : "
                        << Tag::genTag(p, 0, SEND_NODES_TAG) << ")");

      Vector<Real> nodes_to_recv(nodes.storage() + nb_nodes * spatial_dimension,
                                 nb_nodes_to_recv(p) * spatial_dimension);
      comm.receive(nodes_to_recv, p, Tag::genTag(p, 0, SEND_NODES_TAG));
      nb_nodes += nb_nodes_to_recv(p);
    }
#if not defined(AKANTU_NDEBUG)
    else {
      if (p != my_rank) {
        AKANTU_DEBUG_INFO("No nodes to receive from processor " << p);
      }
    }
#endif
  }

  Communicator::waitAll(isend_coordinates_requests);
  Communicator::freeCommunicationRequest(isend_coordinates_requests);

  mesh.sendEvent(new_nodes);
  mesh.sendEvent(new_elements);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template void GridSynchronizer::createGridSynchronizer<IntegrationPoint>(
    const SpatialGrid<IntegrationPoint> & grid);

template void GridSynchronizer::createGridSynchronizer<Element>(
    const SpatialGrid<Element> & grid);

} // namespace akantu
