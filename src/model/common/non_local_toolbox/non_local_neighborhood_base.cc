/**
 * @file   non_local_neighborhood_base.cc
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Mon Sep 21 18:10:49 2015
 *
 * @brief  Implementation of non-local neighborhood base
 *
 * @section LICENSE
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
#include "non_local_neighborhood_base.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
NonLocalNeighborhoodBase::NonLocalNeighborhoodBase(SolidMechanicsModel & model, Real radius)  :
  model(model),
  non_local_radius(radius),
  spatial_grid(NULL), 
  is_creating_grid(false), 
  grid_synchronizer(NULL) {

  AKANTU_DEBUG_IN();


  for(UInt gt = _not_ghost; gt <= _ghost; ++gt) {
    GhostType ghost_type = (GhostType) gt;
    pair_weight[ghost_type] = NULL;
  }

  // if (weight_func_type == "base_wf") 
  //   this->weight_function = new Dummy();
  // case damage_wf:
  //   this->weight_function = new DamagedWeightFunction(); break;
  // case remove_wf:
  //   this->weight_function = new RemoveDamagedWeightFunction(); break;
  // case stress_wf:
  //   this->weight_function = new StressBasedWeightFunction(); break

  this->initNeighborhood();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
NonLocalNeighborhoodBase::~NonLocalNeighborhoodBase() {
  AKANTU_DEBUG_IN();

  for(UInt gt = _not_ghost; gt <= _ghost; ++gt) {
    GhostType ghost_type = (GhostType) gt;
    delete pair_weight[ghost_type];
  }

  delete spatial_grid;
  delete grid_synchronizer;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NonLocalNeighborhoodBase::initNeighborhood() {
  AKANTU_DEBUG_IN();
  //  Material::initMaterial();
  Mesh & mesh = this->model.getMesh();

  ElementTypeMap<UInt> nb_ghost_protected;
  UInt spatial_dimension = this->model.getSpatialDimension();
  Mesh::type_iterator it = mesh.firstType(spatial_dimension, _ghost);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, _ghost);
  for(; it != last_type; ++it)
    nb_ghost_protected(mesh.getNbElement(*it, _ghost), *it, _ghost);

  AKANTU_DEBUG_INFO("Creating the grid");
  this->createGrid();

  /// todo: create pairs, set radius of weight function, intialize weight function, compute weights

  AKANTU_DEBUG_OUT();
}

/* ------------------------------------------------------------------------- */
void NonLocalNeighborhoodBase::createGrid() {
  AKANTU_DEBUG_IN();

  const Real safety_factor = 1.2; // for the cell grid spacing
  Mesh & mesh = this->model.getMesh();
  mesh.computeBoundingBox();

  const Vector<Real> & lower_bounds = mesh.getLocalLowerBounds();
  const Vector<Real> & upper_bounds = mesh.getLocalUpperBounds();
  Vector<Real> center = 0.5 * (upper_bounds + lower_bounds);
  UInt spatial_dimension = this->model.getSpatialDimension();
  Vector<Real> spacing(spatial_dimension, this->non_local_radius * safety_factor);

  spatial_grid = new SpatialGrid<QuadraturePoint>(spatial_dimension, spacing, center);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NonLocalNeighborhoodBase::updatePairList() {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = this->model.getSpatialDimension();
  //// loop over all quads -> all cells
  SpatialGrid<QuadraturePoint>::cells_iterator cell_it = spatial_grid->beginCells();
  SpatialGrid<QuadraturePoint>::cells_iterator cell_end = spatial_grid->endCells();

  Vector<Real> q1_coords(spatial_dimension);
  Vector<Real> q2_coords(spatial_dimension);

  UInt counter = 0;
  for (; cell_it != cell_end; ++cell_it) {
    AKANTU_DEBUG_INFO("Looping on next cell");
    SpatialGrid<QuadraturePoint>::Cell::iterator first_quad =
      spatial_grid->beginCell(*cell_it);
    SpatialGrid<QuadraturePoint>::Cell::iterator last_quad =
      spatial_grid->endCell(*cell_it);
  
    for (;first_quad != last_quad; ++first_quad, ++counter){
      QuadraturePoint q1(*first_quad);
      q1_coords = q1.getPosition();
      AKANTU_DEBUG_INFO("Current quadrature point in this cell: " << q1);
      SpatialGrid<QuadraturePoint>::CellID cell_id = spatial_grid->getCellID(q1_coords);
      /// loop over all the neighbouring cells of the current quad
      SpatialGrid<QuadraturePoint>::neighbor_cells_iterator first_neigh_cell =
  	spatial_grid->beginNeighborCells(cell_id);
      SpatialGrid<QuadraturePoint>::neighbor_cells_iterator last_neigh_cell =
  	spatial_grid->endNeighborCells(cell_id);

      for (; first_neigh_cell != last_neigh_cell; ++first_neigh_cell) {
      	SpatialGrid<QuadraturePoint>::Cell::iterator first_neigh_quad =
      	  spatial_grid->beginCell(*first_neigh_cell);
      	SpatialGrid<QuadraturePoint>::Cell::iterator last_neigh_quad =
      	  spatial_grid->endCell(*first_neigh_cell);

      	// loop over the quadrature point in the current neighboring cell
      	for (;first_neigh_quad != last_neigh_quad; ++first_neigh_quad){
      	  QuadraturePoint q2 = *first_neigh_quad;
      	  q2_coords = q2.getPosition();

      	  Real distance = q1_coords.distance(q2_coords);

      	  if(distance <= this->non_local_radius + Math::getTolerance()  &&
      	     (q2.ghost_type == _ghost ||
      	      (q2.ghost_type == _not_ghost && q1.global_num <= q2.global_num))) { // storing only half lists
      	    pair_list[q2.ghost_type].push_back(std::make_pair(q1, q2));
      	  }
      	}
      }
    }

  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NonLocalNeighborhoodBase::savePairs(const std::string & filename) const {
  std::ofstream pout;

  std::stringstream sstr;

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  Int prank = comm.whoAmI();
  sstr << filename << "." << prank;

  pout.open(sstr.str().c_str());

  for(UInt gt = _not_ghost; gt <= _ghost; ++gt) {
    GhostType ghost_type2 = (GhostType) gt;

    PairList::const_iterator first_pair = pair_list[ghost_type2].begin();
    PairList::const_iterator last_pair  = pair_list[ghost_type2].end();

    for(;first_pair != last_pair; ++first_pair) {

      const QuadraturePoint & lq1 = first_pair->first;
      const QuadraturePoint & lq2 = first_pair->second;
      pout << lq1 << " " << lq2 << " " << std::endl;
    }
  }
}

/* -------------------------------------------------------------------------- */
void NonLocalNeighborhoodBase::saveNeighborCoords(const std::string & filename) const {

  /// this function is not optimazed and only used for tests on small meshes
  /// @todo maybe optimize this function for better performance?

  std::ofstream pout;

  std::stringstream sstr;

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  Int prank = comm.whoAmI();
  sstr << filename << "." << prank;

  pout.open(sstr.str().c_str());

  /// loop over all the quads and write the position of their neighbors
  SpatialGrid<QuadraturePoint>::cells_iterator cell_it = spatial_grid->beginCells();
  SpatialGrid<QuadraturePoint>::cells_iterator cell_end = spatial_grid->endCells();

  for (; cell_it != cell_end; ++cell_it) {
    SpatialGrid<QuadraturePoint>::Cell::iterator first_quad =
      spatial_grid->beginCell(*cell_it);
    SpatialGrid<QuadraturePoint>::Cell::iterator last_quad =
      spatial_grid->endCell(*cell_it);
  
    for (;first_quad != last_quad; ++first_quad){
      QuadraturePoint q1(*first_quad);
      pout << "#neighbors for quad " << q1.global_num << std::endl;
      pout << q1.getPosition() << std::endl;
      for(UInt gt = _not_ghost; gt <= _ghost; ++gt) {
	GhostType ghost_type2 = (GhostType) gt;

	PairList::const_iterator first_pair = pair_list[ghost_type2].begin();
	PairList::const_iterator last_pair  = pair_list[ghost_type2].end();

	for(;first_pair != last_pair; ++first_pair) {
	  if (q1 == first_pair->first && first_pair->second != q1)
	    pout << first_pair->second.getPosition() << std::endl;
	  if (q1 == first_pair->second && first_pair->first != q1)
	    pout << first_pair->first.getPosition() << std::endl;
	}
      }
    }
  }
}


__END_AKANTU__
