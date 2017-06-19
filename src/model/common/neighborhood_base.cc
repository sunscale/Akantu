/**
 * @file   neighborhood_base.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Sat Sep 26 2015
 * @date last modification: Wed Nov 25 2015
 *
 * @brief  Implementation of generic neighborhood base
 *
 * @section LICENSE
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory
 * (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
#include "neighborhood_base.hh"
#include "grid_synchronizer.hh"
#include "model.hh"
/* -------------------------------------------------------------------------- */
#include <fstream>
/* -------------------------------------------------------------------------- */

namespace akantu {
/* -------------------------------------------------------------------------- */
NeighborhoodBase::NeighborhoodBase(Model & model,
                                   const ElementTypeMapReal & quad_coordinates,
                                   const ID & id, const MemoryID & memory_id)
    : Memory(id, memory_id), model(model), neighborhood_radius(0.),
      spatial_grid(NULL), is_creating_grid(false), grid_synchronizer(NULL),
      quad_coordinates(quad_coordinates),
      spatial_dimension(this->model.getMesh().getSpatialDimension()) {

  AKANTU_DEBUG_IN();

  this->registerDataAccessor(*this);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
NeighborhoodBase::~NeighborhoodBase() {
  AKANTU_DEBUG_IN();

  delete spatial_grid;
  delete grid_synchronizer;

  AKANTU_DEBUG_OUT();
}

// /* --------------------------------------------------------------------------
// */
// void NeighborhoodBase::createSynchronizerRegistry(DataAccessor *
// data_accessor){
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
  mesh.computeBoundingBox();

  const Vector<Real> & lower_bounds = mesh.getLocalLowerBounds();
  const Vector<Real> & upper_bounds = mesh.getLocalUpperBounds();
  Vector<Real> center = 0.5 * (upper_bounds + lower_bounds);
  Vector<Real> spacing(spatial_dimension,
                       this->neighborhood_radius * safety_factor);

  spatial_grid =
      new SpatialGrid<IntegrationPoint>(spatial_dimension, spacing, center);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NeighborhoodBase::updatePairList() {
  AKANTU_DEBUG_IN();

  //// loop over all quads -> all cells
  SpatialGrid<IntegrationPoint>::cells_iterator cell_it =
      spatial_grid->beginCells();
  SpatialGrid<IntegrationPoint>::cells_iterator cell_end =
      spatial_grid->endCells();

  Vector<Real> q1_coords(spatial_dimension);
  Vector<Real> q2_coords(spatial_dimension);
  IntegrationPoint q1;
  IntegrationPoint q2;

  UInt counter = 0;
  for (; cell_it != cell_end; ++cell_it) {
    AKANTU_DEBUG_INFO("Looping on next cell");
    SpatialGrid<IntegrationPoint>::Cell::iterator first_quad =
        spatial_grid->beginCell(*cell_it);
    SpatialGrid<IntegrationPoint>::Cell::iterator last_quad =
        spatial_grid->endCell(*cell_it);

    for (; first_quad != last_quad; ++first_quad, ++counter) {
      q1 = *first_quad;
      if (q1.ghost_type == _ghost)
        break;
      Array<Real>::const_vector_iterator coords_type_1_it =
          this->quad_coordinates(q1.type, q1.ghost_type)
              .begin(spatial_dimension);
      q1_coords = coords_type_1_it[q1.global_num];
      AKANTU_DEBUG_INFO("Current quadrature point in this cell: " << q1);
      SpatialGrid<IntegrationPoint>::CellID cell_id =
          spatial_grid->getCellID(q1_coords);
      /// loop over all the neighbouring cells of the current quad
      SpatialGrid<IntegrationPoint>::neighbor_cells_iterator first_neigh_cell =
          spatial_grid->beginNeighborCells(cell_id);
      SpatialGrid<IntegrationPoint>::neighbor_cells_iterator last_neigh_cell =
          spatial_grid->endNeighborCells(cell_id);

      for (; first_neigh_cell != last_neigh_cell; ++first_neigh_cell) {
        SpatialGrid<IntegrationPoint>::Cell::iterator first_neigh_quad =
            spatial_grid->beginCell(*first_neigh_cell);
        SpatialGrid<IntegrationPoint>::Cell::iterator last_neigh_quad =
            spatial_grid->endCell(*first_neigh_cell);

        // loop over the quadrature point in the current neighboring cell
        for (; first_neigh_quad != last_neigh_quad; ++first_neigh_quad) {
          q2 = *first_neigh_quad;
          Array<Real>::const_vector_iterator coords_type_2_it =
              this->quad_coordinates(q2.type, q2.ghost_type)
                  .begin(spatial_dimension);
          q2_coords = coords_type_2_it[q2.global_num];

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
  std::ofstream pout;

  std::stringstream sstr;

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  Int prank = comm.whoAmI();
  sstr << filename << "." << prank;

  pout.open(sstr.str().c_str());

  for (UInt gt = _not_ghost; gt <= _ghost; ++gt) {
    GhostType ghost_type2 = (GhostType)gt;

    PairList::const_iterator first_pair = pair_list[ghost_type2].begin();
    PairList::const_iterator last_pair = pair_list[ghost_type2].end();

    for (; first_pair != last_pair; ++first_pair) {

      const IntegrationPoint & q1 = first_pair->first;
      const IntegrationPoint & q2 = first_pair->second;
      pout << q1 << " " << q2 << " " << std::endl;
    }
  }
}

/* -------------------------------------------------------------------------- */
void NeighborhoodBase::saveNeighborCoords(const std::string & filename) const {

  /// this function is not optimazed and only used for tests on small meshes
  /// @todo maybe optimize this function for better performance?

  Vector<Real> q1_coords(spatial_dimension);
  Vector<Real> q2_coords(spatial_dimension);
  IntegrationPoint q1;
  IntegrationPoint q2;

  std::ofstream pout;

  std::stringstream sstr;

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  Int prank = comm.whoAmI();
  sstr << filename << "." << prank;

  pout.open(sstr.str().c_str());

  /// loop over all the quads and write the position of their neighbors
  SpatialGrid<IntegrationPoint>::cells_iterator cell_it =
      spatial_grid->beginCells();
  SpatialGrid<IntegrationPoint>::cells_iterator cell_end =
      spatial_grid->endCells();

  for (; cell_it != cell_end; ++cell_it) {
    SpatialGrid<IntegrationPoint>::Cell::iterator first_quad =
        spatial_grid->beginCell(*cell_it);
    SpatialGrid<IntegrationPoint>::Cell::iterator last_quad =
        spatial_grid->endCell(*cell_it);

    for (; first_quad != last_quad; ++first_quad) {
      q1 = *first_quad;
      Array<Real>::const_vector_iterator coords_type_1_it =
          this->quad_coordinates(q1.type, q1.ghost_type)
              .begin(spatial_dimension);
      q1_coords = coords_type_1_it[q1.global_num];
      pout << "#neighbors for quad " << q1.global_num << std::endl;
      pout << q1_coords << std::endl;
      for (UInt gt = _not_ghost; gt <= _ghost; ++gt) {
        GhostType ghost_type2 = (GhostType)gt;
        PairList::const_iterator first_pair = pair_list[ghost_type2].begin();
        PairList::const_iterator last_pair = pair_list[ghost_type2].end();
        for (; first_pair != last_pair; ++first_pair) {
          if (q1 == first_pair->first && first_pair->second != q1) {
            q2 = first_pair->second;
            Array<Real>::const_vector_iterator coords_type_2_it =
                this->quad_coordinates(q2.type, q2.ghost_type)
                    .begin(spatial_dimension);
            q2_coords = coords_type_2_it[q2.global_num];
            pout << q2_coords << std::endl;
          }
          if (q1 == first_pair->second && first_pair->first != q1) {
            q2 = first_pair->first;
            Array<Real>::const_vector_iterator coords_type_2_it =
                this->quad_coordinates(q2.type, q2.ghost_type)
                    .begin(spatial_dimension);
            q2_coords = coords_type_2_it[q2.global_num];
            pout << q2_coords << std::endl;
          }
        }
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
void NeighborhoodBase::onElementsRemoved(
    const Array<Element> & element_list,
    const ElementTypeMapArray<UInt> & new_numbering,
    __attribute__((unused)) const RemovedElementsEvent & event) {
  AKANTU_DEBUG_IN();

  FEEngine & fem = this->model.getFEEngine();
  UInt nb_quad = 0;
  // Change the pairs in new global numbering
  for (UInt gt = _not_ghost; gt <= _ghost; ++gt) {
    GhostType ghost_type2 = (GhostType)gt;

    PairList::iterator first_pair = pair_list[ghost_type2].begin();
    PairList::iterator last_pair = pair_list[ghost_type2].end();

    for (; first_pair != last_pair; ++first_pair) {
      IntegrationPoint & q1 = first_pair->first;
      if (new_numbering.exists(q1.type, q1.ghost_type)) {
        UInt q1_new_el = new_numbering(q1.type, q1.ghost_type)(q1.element);
        AKANTU_DEBUG_ASSERT(q1_new_el != UInt(-1), "A local quadrature_point "
                                                   "as been removed instead of "
                                                   "just being renumbered");
        q1.element = q1_new_el;
        nb_quad = fem.getNbIntegrationPoints(q1.type, q1.ghost_type);
        q1.global_num = nb_quad * q1.element + q1.num_point;
      }

      IntegrationPoint & q2 = first_pair->second;
      if (new_numbering.exists(q2.type, q2.ghost_type)) {
        UInt q2_new_el = new_numbering(q2.type, q2.ghost_type)(q2.element);
        AKANTU_DEBUG_ASSERT(q2_new_el != UInt(-1), "A local quadrature_point "
                                                   "as been removed instead of "
                                                   "just being renumbered");
        q2.element = q2_new_el;
        nb_quad = fem.getNbIntegrationPoints(q2.type, q2.ghost_type);
        q2.global_num = nb_quad * q2.element + q2.num_point;
      }
    }
  }
  this->grid_synchronizer->onElementsRemoved(element_list, new_numbering,
                                             event);
  AKANTU_DEBUG_OUT();
}

} // akantu
