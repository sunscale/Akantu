/**
 * @file   material_non_local_inline_impl.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Aug 31 2011
 * @date last modification: Mon Jun 23 2014
 *
 * @brief  Non-local inline implementation
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
__END_AKANTU__

/* -------------------------------------------------------------------------- */
#include "aka_types.hh"
#include "grid_synchronizer.hh"
#include "synchronizer_registry.hh"
#include "integrator.hh"
#include "dumper_paraview.hh"
/* -------------------------------------------------------------------------- */
#include <fstream>
/* -------------------------------------------------------------------------- */

#if defined(AKANTU_DEBUG_TOOLS)
#  include "aka_debug_tools.hh"
#endif


__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt DIM, template <UInt> class WeightFunction>
MaterialNonLocal<DIM, WeightFunction>::MaterialNonLocal(SolidMechanicsModel & model,
							const ID & id)  :
  Material(model, id), weight_func(NULL), spatial_grid(NULL),
  compute_stress_calls(0), is_creating_grid(false), grid_synchronizer(NULL) {
  AKANTU_DEBUG_IN();


  for(UInt gt = _not_ghost; gt <= _ghost; ++gt) {
    GhostType ghost_type = (GhostType) gt;
    pair_weight[ghost_type] = NULL;
  }


  this->is_non_local = true;
  this->weight_func = new WeightFunction<DIM>(*this);

  this->registerSubSection(_st_non_local, "weight_function", *weight_func);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
MaterialNonLocal<spatial_dimension, WeightFunction>::~MaterialNonLocal() {
  AKANTU_DEBUG_IN();

  for(UInt gt = _not_ghost; gt <= _ghost; ++gt) {
    GhostType ghost_type = (GhostType) gt;
    delete pair_weight[ghost_type];
  }

  delete spatial_grid;
  delete weight_func;
  delete grid_synchronizer;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
void MaterialNonLocal<spatial_dimension, WeightFunction>::initMaterial() {
  AKANTU_DEBUG_IN();
  //  Material::initMaterial();
  Mesh & mesh = this->model->getFEEngine().getMesh();

  InternalField<Real> quadrature_points_coordinates("quadrature_points_coordinates_tmp_nl", *this);
  quadrature_points_coordinates.initialize(spatial_dimension);

  ElementTypeMap<UInt> nb_ghost_protected;
  Mesh::type_iterator it = mesh.firstType(spatial_dimension, _ghost);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, _ghost);
  for(; it != last_type; ++it)
    nb_ghost_protected(mesh.getNbElement(*it, _ghost), *it, _ghost);

  AKANTU_DEBUG_INFO("Creating cell list");
  createCellList(quadrature_points_coordinates);

  AKANTU_DEBUG_INFO("Creating pairs");
  updatePairList(quadrature_points_coordinates);

#if not defined(AKANTU_NDEBUG)
  if(AKANTU_DEBUG_TEST(dblDump))
     neighbourhoodStatistics("material_non_local.stats");
#endif

  AKANTU_DEBUG_INFO("Cleaning extra ghosts");
  ///  cleanupExtraGhostElement(nb_ghost_protected);

  AKANTU_DEBUG_INFO("Computing weights");
  weight_func->setRadius(weight_func->getRadius());
  weight_func->init();

  computeWeights(quadrature_points_coordinates);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
void MaterialNonLocal<spatial_dimension, WeightFunction>::cleanupExtraGhostElement(const ElementTypeMap<UInt> & nb_ghost_protected) {
  AKANTU_DEBUG_IN();

  // Create list of element to keep
  std::set<Element> relevant_ghost_element;

  PairList::const_iterator first_pair = pair_list[_ghost].begin();
  PairList::const_iterator last_pair  = pair_list[_ghost].end();
  for(;first_pair != last_pair; ++first_pair) {
    const QuadraturePoint & q2 = first_pair->second;
    relevant_ghost_element.insert(q2);
  }

  // Create list of element to remove and new numbering for element to keep
  Mesh & mesh = this->model->getFEEngine().getMesh();
  std::set<Element> ghost_to_erase;

  Mesh::type_iterator it        = mesh.firstType(spatial_dimension, _ghost);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, _ghost);

  RemovedElementsEvent remove_elem(mesh);


  Element element_local;   // member element corresponds to global element number
  element_local.ghost_type = _ghost;

  for(; it != last_type; ++it) {
    element_local.type = *it;
    UInt nb_ghost_elem_global = mesh.getNbElement(*it, _ghost);
    UInt nb_ghost_elem_protected = 0;
    try {
      nb_ghost_elem_protected = nb_ghost_protected(*it, _ghost);
    } catch (...) {}

    if(!remove_elem.getNewNumbering().exists(*it, _ghost))
      remove_elem.getNewNumbering().alloc(nb_ghost_elem_global, 1, *it, _ghost);
    else remove_elem.getNewNumbering(*it, _ghost).resize(nb_ghost_elem_global);
    Array<UInt> & elem_filter = element_filter(*it, _ghost);
    UInt nb_ghost_elem_local = elem_filter.getSize();
    Array<UInt> & new_numbering = remove_elem.getNewNumbering(*it, _ghost);
    for (UInt g = 0; g < nb_ghost_elem_local; ++g) {
      element_local.element = g;
      Element element_global = this->convertToGlobalElement(element_local);
      if (element_global.element >= nb_ghost_elem_protected &&
	  relevant_ghost_element.find(element_local) == relevant_ghost_element.end()) {
	  // (std::find(relevant_ghost_element.begin(),
	  // 	     relevant_ghost_element.end(),
	  // 	     element_local) == relevant_ghost_element.end())) {
	std::cout<< "element removed" << std::endl;
	remove_elem.getList().push_back(element_global);
	new_numbering(element_global.element) = UInt(-1);
      }
    }

    UInt ng = 0;
    for (UInt g = 0; g < nb_ghost_elem_global; ++g) {
      if (new_numbering(g) != UInt(-1)) {
	new_numbering(g) = ng;
	++ng;
      }
    }
  }

  mesh.sendEvent(remove_elem);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
void MaterialNonLocal<spatial_dimension, WeightFunction>::createCellList(ElementTypeMapArray<Real> & quadrature_points_coordinates) {
  AKANTU_DEBUG_IN();

  const Real safety_factor = 1.2; // for the cell grid spacing
  Mesh & mesh = this->model->getFEEngine().getMesh();
  mesh.computeBoundingBox();

  const Vector<Real> & lower_bounds = mesh.getLocalLowerBounds();
  const Vector<Real> & upper_bounds = mesh.getLocalUpperBounds();
  Vector<Real> center = 0.5 * (upper_bounds + lower_bounds);

  Vector<Real> spacing(spatial_dimension, weight_func->getRadius() * safety_factor);

  spatial_grid = new SpatialGrid<QuadraturePoint>(spatial_dimension, spacing, center);

  this->computeQuadraturePointsCoordinates(quadrature_points_coordinates, _not_ghost);
  this->fillCellList(quadrature_points_coordinates, _not_ghost);
  DumperParaview dumper_ghost("ghosts");
  dumper_ghost.registerMesh(mesh, spatial_dimension, _ghost);
  dumper_ghost.dump();
  is_creating_grid = true;
  std::set<SynchronizationTag> tags;
  tags.insert(_gst_mnl_for_average);
  tags.insert(_gst_mnl_weight);
  tags.insert(_gst_material_id);

  SynchronizerRegistry & synch_registry = this->model->getSynchronizerRegistry();
  std::stringstream sstr; sstr << getID() << ":grid_synchronizer";
  grid_synchronizer = GridSynchronizer::createGridSynchronizer(mesh,
							       *spatial_grid,
							       sstr.str(),
							       &synch_registry,
							       tags);
  is_creating_grid = false;

  this->computeQuadraturePointsCoordinates(quadrature_points_coordinates, _ghost);
  fillCellList(quadrature_points_coordinates, _ghost);
  dumper_ghost.dump();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
void MaterialNonLocal<spatial_dimension, WeightFunction>::fillCellList(const ElementTypeMapArray<Real> & quadrature_points_coordinates,
								       const GhostType & ghost_type) {
  QuadraturePoint q;
  q.ghost_type = ghost_type;

  Mesh::type_iterator it        = this->element_filter.firstType(spatial_dimension, ghost_type);
  Mesh::type_iterator last_type = this->element_filter.lastType (spatial_dimension, ghost_type);
  for(; it != last_type; ++it) {
    Array<UInt> & elem_filter = element_filter(*it, ghost_type);
    UInt nb_element = elem_filter.getSize();
    UInt nb_quad    = this->model->getFEEngine().getNbQuadraturePoints(*it, ghost_type);

    const Array<Real> & quads = quadrature_points_coordinates(*it, ghost_type);
    q.type = *it;

    Array<Real>::const_vector_iterator quad = quads.begin(spatial_dimension);
    UInt * elem = elem_filter.storage();

    for (UInt e = 0; e < nb_element; ++e) {
      q.element = *elem;
      for (UInt nq = 0; nq < nb_quad; ++nq) {
	q.num_point = nq;
	q.global_num = q.element * nb_quad + nq;
	//q.setPosition(*quad);
	spatial_grid->insert(q, *quad);
	++quad;
      }
      ++elem;
    }
  }
}


/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
void MaterialNonLocal<spatial_dimension, WeightFunction>::updatePairList(const ElementTypeMapArray<Real> & quadrature_points_coordinates) {
  AKANTU_DEBUG_IN();

  GhostType ghost_type = _not_ghost;
  QuadraturePoint q1;
  q1.ghost_type = ghost_type;

  // generate the pair of neighbor depending of the cell_list
  Mesh::type_iterator it        = this->element_filter.firstType(spatial_dimension, ghost_type);
  Mesh::type_iterator last_type = this->element_filter.lastType(spatial_dimension, ghost_type);
  for(; it != last_type; ++it) {
    // Preparing datas
    const Array<Real> & quads = quadrature_points_coordinates(*it, ghost_type);
    Array<Real>::const_vector_iterator q1_coord_it = quads.begin(spatial_dimension);
    Array<Real>::const_vector_iterator last_quad   = quads.end(spatial_dimension);

    q1.type = *it;
    q1.global_num = 0;
    // loop over quad points
    for(;q1_coord_it != last_quad; ++q1_coord_it, ++(q1.global_num)) {
      UInt nb_quad1 =
	this->model->getFEEngine().getNbQuadraturePoints(q1.type,
							 q1.ghost_type);
      q1.element   = q1.global_num / nb_quad1;
      q1.num_point = q1.global_num % nb_quad1;
      const Vector<Real> & q1_coord = *q1_coord_it;

      SpatialGrid<QuadraturePoint>::CellID cell_id = spatial_grid->getCellID(q1_coord);
      SpatialGrid<QuadraturePoint>::neighbor_cells_iterator first_neigh_cell =
	spatial_grid->beginNeighborCells(cell_id);
      SpatialGrid<QuadraturePoint>::neighbor_cells_iterator last_neigh_cell =
	spatial_grid->endNeighborCells(cell_id);

      // loop over neighbors cells of the one containing the current quadrature
      // point
      for (; first_neigh_cell != last_neigh_cell; ++first_neigh_cell) {
	SpatialGrid<QuadraturePoint>::Cell::iterator first_neigh_quad =
	  spatial_grid->beginCell(*first_neigh_cell);
	SpatialGrid<QuadraturePoint>::Cell::iterator last_neigh_quad =
	  spatial_grid->endCell(*first_neigh_cell);

	// loop over the quadrature point in the current cell of the cell list
	for (;first_neigh_quad != last_neigh_quad; ++first_neigh_quad){
	  QuadraturePoint q2 = this->convertToLocalPoint(*first_neigh_quad);

	  Array<Real>::const_vector_iterator q2_coord_it =
	    quadrature_points_coordinates(q2.type,
					  q2.ghost_type).begin(spatial_dimension);

	  const Vector<Real> & q2_coord = q2_coord_it[q2.global_num];

	  Real distance = q1_coord.distance(q2_coord);

	  if(distance <= weight_func->getRadius() &&
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
template<UInt spatial_dimension, template <UInt> class WeightFunction>
void MaterialNonLocal<spatial_dimension, WeightFunction>::computeWeights(const ElementTypeMapArray<Real> & quadrature_points_coordinates) {
  AKANTU_DEBUG_IN();

  InternalField<Real> quadrature_points_volumes("quadrature_points_volumes", *this);
  quadrature_points_volumes.initialize(1);

  const FEEngine & fem = this->model->getFEEngine();

  weight_func->updateInternals(quadrature_points_volumes);

  for(UInt gt = _not_ghost; gt <= _ghost; ++gt) {
    GhostType ghost_type2 = (GhostType) gt;

    if(!(pair_weight[ghost_type2])) {
      std::string ghost_id = "";
      if (ghost_type2 == _ghost) ghost_id = ":ghost";
      std::stringstream sstr; sstr << getID() << ":pair_weight:" << ghost_id;
      pair_weight[ghost_type2] = new Array<Real>(0, 2, sstr.str());
    }

    pair_weight[ghost_type2]->resize(pair_list[ghost_type2].size());
    pair_weight[ghost_type2]->clear();

    PairList::const_iterator first_pair = pair_list[ghost_type2].begin();
    PairList::const_iterator last_pair  = pair_list[ghost_type2].end();

    Array<Real>::vector_iterator weight_it = pair_weight[ghost_type2]->begin(2);

    // Compute the weights
    for(;first_pair != last_pair; ++first_pair, ++weight_it) {
      Vector<Real> & weight = *weight_it;

      const QuadraturePoint & lq1 = first_pair->first;
      const QuadraturePoint & lq2 = first_pair->second;

      QuadraturePoint gq1 = this->convertToGlobalPoint(lq1);
      QuadraturePoint gq2 = this->convertToGlobalPoint(lq2);

      //   const Real q2_wJ = fem.getIntegratorInterface().getJacobians(gq2.type, gq2.ghost_type)(gq2.global_num);

      Array<Real>::const_vector_iterator quad_coords_1 =
	quadrature_points_coordinates(lq1.type, lq1.ghost_type).begin(spatial_dimension);
      Array<Real>::const_vector_iterator quad_coords_2 =
	quadrature_points_coordinates(lq2.type, lq2.ghost_type).begin(spatial_dimension);

      Array<Real> & quad_volumes_1 = quadrature_points_volumes(lq1.type, lq1.ghost_type);
      const Array<Real> & jacobians_2 =
	fem.getIntegratorInterface().getJacobians(gq2.type, gq2.ghost_type);
      const Real & q2_wJ = jacobians_2(gq2.global_num);
      // Real & q1_volume = quad_volumes(lq1.global_num);

      const Vector<Real> & q1_coord = quad_coords_1[lq1.global_num];
      // quadrature_points_coordinates(lq1.type, lq1.ghost_type).begin(spatial_dimension)[lq1.global_num];
      const Vector<Real> & q2_coord = quad_coords_2[lq2.global_num];
      // quadrature_points_coordinates(lq1.type, lq1.ghost_type).begin(spatial_dimension)[lq2.global_num];

      this->weight_func->selectType(lq1.type, lq1.ghost_type, lq2.type, lq2.ghost_type);

      // Weight function
      Real r = q1_coord.distance(q2_coord);
      Real w1 = this->weight_func->operator()(r, lq1, lq2);
      weight(0) = q2_wJ * w1;
      //     q1_volume += weight(0);
      quad_volumes_1(lq1.global_num) += weight(0);

      if(lq2.ghost_type != _ghost && lq1.global_num != lq2.global_num) {
	const Array<Real> & jacobians_1 =
	  fem.getIntegratorInterface().getJacobians(gq1.type, gq1.ghost_type);
	Array<Real> & quad_volumes_2 =
	  quadrature_points_volumes(lq2.type, lq2.ghost_type);

	const Real & q1_wJ = jacobians_1(gq1.global_num);
	//Real & q2_volume = quad_volumes_2(lq2.global_num);

	Real w2 = this->weight_func->operator()(r, lq2, lq1);
	weight(1) = q1_wJ * w2;

	quad_volumes_2(lq2.global_num) += weight(1);
	//q2_volume += weight(1);
      } else
	weight(1) = 0.;
    }
  }

  //normalize the weights
  for(UInt gt = _not_ghost; gt <= _ghost; ++gt) {
    GhostType ghost_type2 = (GhostType) gt;

    PairList::const_iterator first_pair = pair_list[ghost_type2].begin();
    PairList::const_iterator last_pair  = pair_list[ghost_type2].end();

    Array<Real>::vector_iterator weight_it = pair_weight[ghost_type2]->begin(2);

    // Compute the weights
    for(;first_pair != last_pair; ++first_pair, ++weight_it) {
      Vector<Real> & weight = *weight_it;

      const QuadraturePoint & lq1 = first_pair->first;
      const QuadraturePoint & lq2 = first_pair->second;

      Array<Real> & quad_volumes_1 = quadrature_points_volumes(lq1.type, lq1.ghost_type);
      Array<Real> & quad_volumes_2 = quadrature_points_volumes(lq2.type, lq2.ghost_type);

      Real q1_volume = quad_volumes_1(lq1.global_num);

      weight(0) *= 1. / q1_volume;
      if(ghost_type2 != _ghost) {
	Real q2_volume = quad_volumes_2(lq2.global_num);
	weight(1) *= 1. / q2_volume;
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
template<typename T>
void MaterialNonLocal<spatial_dimension, WeightFunction>::weightedAvergageOnNeighbours(const InternalField<T> & to_accumulate,
										       InternalField<T> & accumulated,
										       UInt nb_degree_of_freedom,
										       GhostType ghost_type2) const {
  AKANTU_DEBUG_IN();

  if(ghost_type2 == _not_ghost)  {
    accumulated.reset();
  }

  PairList::const_iterator first_pair = pair_list[ghost_type2].begin();
  PairList::const_iterator last_pair  = pair_list[ghost_type2].end();

  Array<Real>::vector_iterator weight_it = pair_weight[ghost_type2]->begin(2);

  // Compute the weights
  for(;first_pair != last_pair; ++first_pair, ++weight_it) {
    Vector<Real> & weight = *weight_it;

    const QuadraturePoint & lq1 = first_pair->first;
    const QuadraturePoint & lq2 = first_pair->second;

    const Array<T> & to_acc_1 = to_accumulate(lq1.type, lq1.ghost_type);
    Array<T> & acc_1 = accumulated(lq1.type, lq1.ghost_type);
    const Array<T> & to_acc_2 = to_accumulate(lq2.type, lq2.ghost_type);
    Array<T> & acc_2 = accumulated(lq2.type, lq2.ghost_type);

    // const Vector<T> & q2_to_acc = to_acc_2[lq2.global_num];
    // Vector<T> & q1_acc = acc_1[lq1.global_num];

    //q1_acc += weight(0) * q2_to_acc;
    for(UInt d = 0; d < nb_degree_of_freedom; ++d) {
      acc_1(lq1.global_num, d) += weight(0) * to_acc_2(lq2.global_num, d);
    }

    if(ghost_type2 != _ghost) {
      for(UInt d = 0; d < nb_degree_of_freedom; ++d) {
      	acc_2(lq2.global_num, d) += weight(1) * to_acc_1(lq1.global_num, d);
      }
    }

    // if(ghost_type2 != _ghost) {
    //   const Vector<T> & q1_to_acc = to_acc_1[lq1.global_num];
    //   Vector<T> & q2_acc = acc_2[lq2.global_num];

    //   q2_acc += weight(1) * q1_to_acc;
    // }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
void MaterialNonLocal<spatial_dimension, WeightFunction>::updateResidual(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  // Update the weights for the non local variable averaging
  if(ghost_type == _not_ghost &&
     this->weight_func->getUpdateRate() &&
     (this->compute_stress_calls % this->weight_func->getUpdateRate() == 0)) {
    ElementTypeMapArray<Real> quadrature_points_coordinates("quadrature_points_coordinates", getID());
    Mesh & mesh = this->model->getFEEngine().getMesh();
    mesh.initElementTypeMapArray(quadrature_points_coordinates, spatial_dimension, spatial_dimension);
    computeQuadraturePointsCoordinates(quadrature_points_coordinates, _not_ghost);
    computeQuadraturePointsCoordinates(quadrature_points_coordinates, _ghost);
    computeWeights(quadrature_points_coordinates);
  }
  if(ghost_type == _not_ghost) ++this->compute_stress_calls;

  computeAllStresses(ghost_type);

  computeNonLocalStresses(ghost_type);
  assembleResidual(ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
void MaterialNonLocal<spatial_dimension, WeightFunction>::computeAllNonLocalStresses(GhostType ghost_type) {
  // Update the weights for the non local variable averaging
  if(ghost_type == _not_ghost) {
    if(this->weight_func->getUpdateRate() &&
       (this->compute_stress_calls % this->weight_func->getUpdateRate() == 0)) {
      this->model->getSynchronizerRegistry().asynchronousSynchronize(_gst_mnl_weight);

      ElementTypeMapArray<Real> quadrature_points_coordinates("quadrature_points_coordinates", getID());
      Mesh & mesh = this->model->getFEEngine().getMesh();
      mesh.initElementTypeMapArray(quadrature_points_coordinates, spatial_dimension, spatial_dimension);
      computeQuadraturePointsCoordinates(quadrature_points_coordinates, _not_ghost);
      computeQuadraturePointsCoordinates(quadrature_points_coordinates, _ghost);

      this->model->getSynchronizerRegistry().waitEndSynchronize(_gst_mnl_weight);

      computeWeights(quadrature_points_coordinates);
    }

    typename std::map<ID, NonLocalVariable>::iterator it = non_local_variables.begin();
    typename std::map<ID, NonLocalVariable>::iterator end = non_local_variables.end();
    for(;it != end; ++it) {
      NonLocalVariable & non_local_variable = it->second;

      non_local_variable.non_local->resize();
      this->weightedAvergageOnNeighbours(*non_local_variable.local, *non_local_variable.non_local,
					 non_local_variable.nb_component, _not_ghost);
    }

    ++this->compute_stress_calls;
  } else {
    typename std::map<ID, NonLocalVariable>::iterator it = non_local_variables.begin();
    typename std::map<ID, NonLocalVariable>::iterator end = non_local_variables.end();
    for(;it != end; ++it) {
      NonLocalVariable & non_local_variable = it->second;
      this->weightedAvergageOnNeighbours(*non_local_variable.local, *non_local_variable.non_local,
					 non_local_variable.nb_component, _ghost);
    }
    computeNonLocalStresses(_not_ghost);
  }
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
void MaterialNonLocal<spatial_dimension, WeightFunction>::savePairs(const std::string & filename) const {
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

    Array<Real>::vector_iterator weight_it = pair_weight[ghost_type2]->begin(2);

    for(;first_pair != last_pair; ++first_pair, ++weight_it) {
      Vector<Real> & weight = *weight_it;

      const QuadraturePoint & lq1 = first_pair->first;
      const QuadraturePoint & lq2 = first_pair->second;
      pout << lq1 << " " << lq2 << " " << weight << std::endl;
    }
  }
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
void MaterialNonLocal<spatial_dimension, WeightFunction>::neighbourhoodStatistics(const std::string & filename) const {
  //   std::ofstream pout;
  // pout.open(filename.c_str());

  // const Mesh & mesh = this->model->getFEEngine().getMesh();

  // GhostType ghost_type1;
  // ghost_type1 = _not_ghost;

  // StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  // Int prank = comm.whoAmI();

  // InternalField<UInt> nb_neighbors("nb_neighbours", *const_cast<MaterialNonLocal *>(this));
  // nb_neighbors.initialize(1);

  // for(UInt gt = _not_ghost; gt <= _ghost; ++gt) {
  //   GhostType ghost_type2 = (GhostType) gt;
  //   UInt existing_pairs_num = gt - _not_ghost;

  //   pair_type::const_iterator first_pair_types = existing_pairs[existing_pairs_num].begin();
  //   pair_type::const_iterator last_pair_types = existing_pairs[existing_pairs_num].end();

  //   for (; first_pair_types != last_pair_types; ++first_pair_types) {
  //     const Array<UInt> & pairs =
  //	pair_list(first_pair_types->first, ghost_type1)(first_pair_types->second, ghost_type2);
  //     if(prank == 0) {
  //	pout << ghost_type2 << " ";
  //	pout << "Types : " << first_pair_types->first << " " << first_pair_types->second << std::endl;
  //     }
  //     Array<UInt>::const_iterator< Vector<UInt> > first_pair = pairs.begin(2);
  //     Array<UInt>::const_iterator< Vector<UInt> > last_pair  = pairs.end(2);
  //     Array<UInt> & nb_neigh_1 = nb_neighbors(first_pair_types->first, ghost_type1);
  //     Array<UInt> & nb_neigh_2 = nb_neighbors(first_pair_types->second, ghost_type2);
  //     for(;first_pair != last_pair; ++first_pair) {
  //	UInt q1 = (*first_pair)(0);
  //	UInt q2 = (*first_pair)(1);
  //	++(nb_neigh_1(q1));
  //	if(q1 != q2) ++(nb_neigh_2(q2));
  //     }
  //   }

  //   Mesh::type_iterator it        = mesh.firstType(spatial_dimension, ghost_type1);
  //   Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, ghost_type1);
  //   UInt nb_quads = 0;
  //   Real sum_nb_neig = 0;
  //   UInt max_nb_neig = 0;
  //   UInt min_nb_neig = std::numeric_limits<UInt>::max();
  //   for(; it != last_type; ++it) {
  //     Array<UInt> & nb_neighor = nb_neighbors(*it, ghost_type1);
  //     Array<UInt>::iterator<UInt> nb_neigh = nb_neighor.begin();
  //     Array<UInt>::iterator<UInt> end_neigh  = nb_neighor.end();

  //     for (; nb_neigh != end_neigh; ++nb_neigh, ++nb_quads) {
  //	UInt nb = *nb_neigh;
  //	sum_nb_neig += nb;
  //	max_nb_neig = std::max(max_nb_neig, nb);
  //	min_nb_neig = std::min(min_nb_neig, nb);
  //     }
  //   }


  //   comm.allReduce(&nb_quads,    1, _so_sum);
  //   comm.allReduce(&sum_nb_neig, 1, _so_sum);
  //   comm.allReduce(&max_nb_neig, 1, _so_max);
  //   comm.allReduce(&min_nb_neig, 1, _so_min);

  //   if(prank == 0) {
  //     pout << ghost_type2 << " ";
  //     pout << "Nb quadrature points: " << nb_quads << std::endl;

  //     Real mean_nb_neig = sum_nb_neig / Real(nb_quads);
  //     pout << ghost_type2 << " ";
  //     pout << "Average nb neighbors: " << mean_nb_neig << "(" << sum_nb_neig << ")" << std::endl;

  //     pout << ghost_type2 << " ";
  //     pout << "Max nb neighbors:     " << max_nb_neig << std::endl;

  //     pout << ghost_type2 << " ";
  //     pout << "Min nb neighbors:     " << min_nb_neig << std::endl;
  //   }
  // }
  // pout.close();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
inline UInt MaterialNonLocal<spatial_dimension, WeightFunction>::getNbDataForElements(const Array<Element> & elements,
										      SynchronizationTag tag) const {
  UInt nb_quadrature_points = this->getModel().getNbQuadraturePoints(elements);
  UInt size = 0;

  if(tag == _gst_mnl_for_average) {
    typename std::map<ID, NonLocalVariable>::const_iterator it = non_local_variables.begin();
    typename std::map<ID, NonLocalVariable>::const_iterator end = non_local_variables.end();

    for(;it != end; ++it) {
      const NonLocalVariable & non_local_variable = it->second;
      size += non_local_variable.nb_component * sizeof(Real) * nb_quadrature_points;
    }
  }

  size += weight_func->getNbDataForElements(elements, tag);

  return size;
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
inline void MaterialNonLocal<spatial_dimension, WeightFunction>::packElementData(CommunicationBuffer & buffer,
										 const Array<Element> & elements,
										 SynchronizationTag tag) const {
  if(tag == _gst_mnl_for_average) {
    typename std::map<ID, NonLocalVariable>::const_iterator it = non_local_variables.begin();
    typename std::map<ID, NonLocalVariable>::const_iterator end = non_local_variables.end();

    for(;it != end; ++it) {
      const NonLocalVariable & non_local_variable = it->second;
      this->packElementDataHelper(*non_local_variable.local,
				  buffer, elements);
    }
  }

  weight_func->packElementData(buffer, elements, tag);
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
inline void MaterialNonLocal<spatial_dimension, WeightFunction>::unpackElementData(CommunicationBuffer & buffer,
										   const Array<Element> & elements,
										   SynchronizationTag tag) {
  if(tag == _gst_mnl_for_average) {
    typename std::map<ID, NonLocalVariable>::iterator it = non_local_variables.begin();
    typename std::map<ID, NonLocalVariable>::iterator end = non_local_variables.end();

    for(;it != end; ++it) {
      NonLocalVariable & non_local_variable = it->second;
      this->unpackElementDataHelper(*non_local_variable.local,
				    buffer, elements);
    }
  }

  weight_func->unpackElementData(buffer, elements, tag);
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
inline void MaterialNonLocal<spatial_dimension, WeightFunction>::onElementsAdded(const Array<Element> & element_list) {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ERROR("This is a case not taken into account!!!");
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
inline void MaterialNonLocal<spatial_dimension, WeightFunction>::onElementsRemoved(const Array<Element> & element_list,
										   const ElementTypeMapArray<UInt> & new_numbering,
										   __attribute__((unused)) const RemovedElementsEvent & event) {
  AKANTU_DEBUG_IN();

  // Change the pairs in new global numbering
  for(UInt gt = _not_ghost; gt <= _ghost; ++gt) {
    GhostType ghost_type2 = (GhostType) gt;

    PairList::iterator first_pair = pair_list[ghost_type2].begin();
    PairList::iterator last_pair  = pair_list[ghost_type2].end();

    //   Array<Real>::vector_iterator weight_it = pair_weight[ghost_type2]->begin(2);

    for(;first_pair != last_pair; ++first_pair) {
      QuadraturePoint & q1 = first_pair->first;
      QuadraturePoint gq1  = this->convertToGlobalPoint(q1);
      q1 = gq1;

      if(new_numbering.exists(q1.type, q1.ghost_type)) {
	UInt q1_new_el = new_numbering(q1.type, q1.ghost_type)(gq1.element);
	AKANTU_DEBUG_ASSERT(q1_new_el != UInt(-1), "A local quadrature_point as been removed instead of just being renumbered");
	q1.element = q1_new_el;
      }


      QuadraturePoint & q2 = first_pair->second;
      QuadraturePoint gq2  = this->convertToGlobalPoint(q2);
      q2 = gq2;

      if(new_numbering.exists(q2.type, q2.ghost_type)) {
	UInt q2_new_el = new_numbering(q2.type, q2.ghost_type)(gq2.element);
	AKANTU_DEBUG_ASSERT(q2_new_el != UInt(-1), "A local quadrature_point as been removed instead of just being renumbered");
	q2.element = q2_new_el;
      }
    }
  }

  // Change the material numbering
  Material::onElementsRemoved(element_list, new_numbering, event);

  // Change back the pairs to the new material numbering
  for(UInt gt = _not_ghost; gt <= _ghost; ++gt) {
    GhostType ghost_type2 = (GhostType) gt;

    PairList::iterator first_pair = pair_list[ghost_type2].begin();
    PairList::iterator last_pair  = pair_list[ghost_type2].end();

    //   Array<Real>::vector_iterator weight_it = pair_weight[ghost_type2]->begin(2);

    for(;first_pair != last_pair; ++first_pair) {
      first_pair->first  = this->convertToLocalPoint(first_pair->first );
      first_pair->second = this->convertToLocalPoint(first_pair->second);
    }
  }

  AKANTU_DEBUG_OUT();
}
