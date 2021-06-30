/**
 * @file   neighborhood_base.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Sat Sep 26 2015
 * @date last modification: Wed Jan 31 2018
 *
 * @brief  Implementation of generic neighborhood base
 *
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "neighborhood_base.hh"
#include "grid_synchronizer.hh"
#include "mesh_accessor.hh"
#include "model.hh"
/* -------------------------------------------------------------------------- */
#include <fstream>
/* -------------------------------------------------------------------------- */

namespace akantu {
/* -------------------------------------------------------------------------- */
NeighborhoodBase::NeighborhoodBase(Model & model,
                                   const ElementTypeMapReal & quad_coordinates,
                                   const ID & id)
    : id(id), model(model), quad_coordinates(quad_coordinates),
      spatial_dimension(this->model.getMesh().getSpatialDimension()) {

  AKANTU_DEBUG_IN();

  this->registerDataAccessor(*this);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
NeighborhoodBase::~NeighborhoodBase() = default;

/* -------------------------------------------------------------------------- */
// void NeighborhoodBase::createSynchronizerRegistry(
//     DataAccessor<Element> * data_accessor) {
//   this->synch_registry = new SynchronizerRegistry(*data_accessor);
// }

/* -------------------------------------------------------------------------- */
void NeighborhoodBase::initNeighborhood() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Creating the grid");
  this->createGrid();

  AKANTU_DEBUG_OUT();
}

/* ------------------------------------------------------------------------- */
void NeighborhoodBase::createGrid() {
  AKANTU_DEBUG_IN();

  const Real safety_factor = 1.2; // for the cell grid spacing
  Mesh & mesh = this->model.getMesh();

  const auto & lower_bounds = mesh.getLocalLowerBounds();
  const auto & upper_bounds = mesh.getLocalUpperBounds();
  Vector<Real> center = 0.5 * (upper_bounds + lower_bounds);
  Vector<Real> spacing(spatial_dimension,
                       this->neighborhood_radius * safety_factor);

  spatial_grid = std::make_unique<SpatialGrid<IntegrationPoint>>(
      spatial_dimension, spacing, center);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NeighborhoodBase::updatePairList() {
  AKANTU_DEBUG_IN();

  //// loop over all quads -> all cells
  for (auto && cell_id : *spatial_grid) {
    AKANTU_DEBUG_INFO("Looping on next cell");

    for (auto && q1 : spatial_grid->getCell(cell_id)) {
      if (q1.ghost_type == _ghost) {
        break;
      }
      auto coords_type_1_it = this->quad_coordinates(q1.type, q1.ghost_type)
                                  .begin(spatial_dimension);
      auto q1_coords = Vector<Real>(coords_type_1_it[q1.global_num]);

      AKANTU_DEBUG_INFO("Current quadrature point in this cell: " << q1);
      auto cell_id = spatial_grid->getCellID(q1_coords);

      /// loop over all the neighboring cells of the current quad
      for (auto && neighbor_cell : cell_id.neighbors()) {
        // loop over the quadrature point in the current neighboring cell
        for (auto && q2 : spatial_grid->getCell(neighbor_cell)) {
          auto coords_type_2_it = this->quad_coordinates(q2.type, q2.ghost_type)
                                      .begin(spatial_dimension);
          auto q2_coords = Vector<Real>(coords_type_2_it[q2.global_num]);

          Real distance = q1_coords.distance(q2_coords);

          if (distance <= this->neighborhood_radius + Math::getTolerance() &&
              (q2.ghost_type == _ghost ||
               (q2.ghost_type == _not_ghost &&
                q1.global_num <= q2.global_num))) { // storing only half lists
            pair_list[q2.ghost_type].push_back(std::make_pair(q1, q2));
          }
        }
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NeighborhoodBase::savePairs(const std::string & filename) const {
  std::stringstream sstr;

  const Communicator & comm = model.getMesh().getCommunicator();
  Int prank = comm.whoAmI();
  sstr << filename << "." << prank;

  std::ofstream pout;
  pout.open(sstr.str().c_str());

  for (auto && ghost_type : ghost_types) {
    for (const auto & pair : pair_list[ghost_type]) {
      const auto & q1 = pair.first;
      const auto & q2 = pair.second;
      pout << q1 << " " << q2 << " " << std::endl;
    }
  }

  pout.close();

  if (comm.getNbProc() != 1) {
    return;
  }

  Mesh mesh_out(spatial_dimension);
  MeshAccessor mesh_accessor(mesh_out);
  auto & connectivity = mesh_accessor.getConnectivity(_segment_2);
  auto & tag = mesh_accessor.getData<UInt>("tag_1", _segment_2);
  auto & nodes = mesh_accessor.getNodes();

  std::map<IntegrationPoint, UInt> quad_to_nodes;
  UInt node = 0;

  IntegrationPoint q1;
  IntegrationPoint q2;
  bool inserted;
  for (auto && ghost_type : ghost_types) {
    for (const auto & pair : pair_list[ghost_type]) {
      std::tie(q1, q2) = pair;

      auto add_node = [&](auto && q) {
        std::tie(std::ignore, inserted) =
            quad_to_nodes.insert(std::make_pair(q, node));

        if (not inserted) {
          return;
        }

        auto coords_it = this->quad_coordinates(q.type, q.ghost_type)
                             .begin(spatial_dimension);
        auto && coords = Vector<Real>(coords_it[q.global_num]);
        nodes.push_back(coords);
        ++node;
      };

      add_node(q1);
      add_node(q2);
    }
  }

  for (auto && ghost_type : ghost_types) {
    for (const auto & pair : pair_list[ghost_type]) {
      std::tie(q1, q2) = pair;

      UInt node1 = quad_to_nodes[q1];
      UInt node2 = quad_to_nodes[q2];

      connectivity.push_back(Vector<UInt>{node1, node2});
      tag.push_back(node1 + 1);
      if (node1 != node2) {
        connectivity.push_back(Vector<UInt>{node2, node1});
        tag.push_back(node2 + 1);
      }
    }
  }

  mesh_out.write(filename + ".msh");
}

/* -------------------------------------------------------------------------- */
void NeighborhoodBase::saveNeighborCoords(const std::string & filename) const {
  // this function is not optimized and only used for tests on small meshes
  // @todo maybe optimize this function for better performance?
  IntegrationPoint q2;

  std::stringstream sstr;

  const Communicator & comm = model.getMesh().getCommunicator();
  Int prank = comm.whoAmI();
  sstr << filename << "." << prank;

  std::ofstream pout;
  pout.open(sstr.str().c_str());

  /// loop over all the quads and write the position of their neighbors
  for (auto && cell_id : *spatial_grid) {
    for (auto && q1 : spatial_grid->getCell(cell_id)) {
      auto coords_type_1_it = this->quad_coordinates(q1.type, q1.ghost_type)
                                  .begin(spatial_dimension);
      auto && q1_coords = Vector<Real>(coords_type_1_it[q1.global_num]);

      pout << "#neighbors for quad " << q1.global_num << std::endl;
      pout << q1_coords << std::endl;

      for (auto && ghost_type2 : ghost_types) {
        for (auto && pair : pair_list[ghost_type2]) {
          if (q1 == pair.first && pair.second != q1) {
            q2 = pair.second;
          } else if (q1 == pair.second && pair.first != q1) {
            q2 = pair.first;
          } else {
            continue;
          }

          auto coords_type_2_it = this->quad_coordinates(q2.type, q2.ghost_type)
                                      .begin(spatial_dimension);
          auto && q2_coords = Vector<Real>(coords_type_2_it[q2.global_num]);
          pout << q2_coords << std::endl;
        }
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
void NeighborhoodBase::onElementsRemoved(
    const Array<Element> & element_list,
    const ElementTypeMapArray<UInt> & new_numbering,
    const RemovedElementsEvent & event) {
  AKANTU_DEBUG_IN();

  FEEngine & fem = this->model.getFEEngine();
  UInt nb_quad = 0;
  auto cleanPoint = [&](auto && q) {
    if (new_numbering.exists(q.type, q.ghost_type)) {
      UInt q_new_el = new_numbering(q.type, q.ghost_type)(q.element);
      AKANTU_DEBUG_ASSERT(q_new_el != UInt(-1),
                          "A local quadrature_point "
                              << q
                              << " as been removed instead of "
                                 "just being renumbered: "
                              << id);
      q.element = q_new_el;
      nb_quad = fem.getNbIntegrationPoints(q.type, q.ghost_type);
      q.global_num = nb_quad * q.element + q.num_point;
    }
  };

  // Change the pairs in new global numbering
  for (auto ghost_type : ghost_types) {
    auto & pair_list = this->pair_list.at(ghost_type);
    for (auto && pair : pair_list) {
      if (pair.first.ghost_type == _ghost) {
        cleanPoint(pair.first);
      }
      if (pair.second.ghost_type == _ghost) {
        cleanPoint(pair.second);
      }
    }
  }

  this->grid_synchronizer->onElementsRemoved(element_list, new_numbering,
                                             event);
  AKANTU_DEBUG_OUT();
}

} // namespace akantu
