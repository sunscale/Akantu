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

  this->is_non_local = true;
  this->weight_func = new WeightFunction<DIM>(*this);

  this->registerSubSection(_st_non_local, "weight_function", *weight_func);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
MaterialNonLocal<spatial_dimension, WeightFunction>::~MaterialNonLocal() {
  AKANTU_DEBUG_IN();

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
  cleanupExtraGhostElement(nb_ghost_protected);

  AKANTU_DEBUG_INFO("Computing weights");
  weight_func->setRadius(weight_func->getRadius());
  weight_func->init();

  computeWeights(quadrature_points_coordinates);

  // >>>>>> DEBUG CODE >>>>>> //
#if defined(AKANTU_DEBUG_TOOLS) && defined(AKANTU_CORE_CXX11)
  debug::element_manager.print
    (debug::_dm_material_non_local,
     [this, &quadrature_points_coordinates](const Element & el)->std::string {
      std::stringstream out;
      FEEngine & fem = this->model->getFEEngine();

      Mesh mesh(spatial_dimension, getID() + "mesh_tmp");
      mesh.addConnectivityType(_segment_2);
      Array<UInt> & connectivity = const_cast<Array<UInt> &>(mesh.getConnectivity(_segment_2));
      Array<Real> & nodes = const_cast<Array<Real> &>(mesh.getNodes());
      Array<bool> g_notg(0,1);
      Array<Real> dist(0,1);
      Array<Real> weight(0,1);
      Array<Real> damage(0,1);
      Array<Real> jac(0,1);

      UInt nb_quad_per_elem =
        this->model->getFEEngine().getNbQuadraturePoints(el.type,
                                                    el.ghost_type);

      Array<Real>::const_vector_iterator quad_it = quadrature_points_coordinates(el.type, el.ghost_type).begin(spatial_dimension);
      std::map<Vector<Real>, UInt> numbering;
      UInt counter = 0;
      const Vector<Real> quad = quad_it[el.element * nb_quad_per_elem];
      numbering[el.element * nb_quad_per_elem] = counter++;
      nodes.push_back(quad);
      g_notg.push_back(el.ghost_type == _ghost);
      dist.push_back(quad.distance(quad));
      weight.push_back(0.);
      damage.push_back((getArray("damage", el.type, el.ghost_type)(el.element * fem.getNbQuadraturePoints(el.type, el.ghost_type))));
      jac.push_back(fem.getIntegratorInterface().getJacobians(el.type, el.ghost_type)(el.element * fem.getNbQuadraturePoints(el.type, el.ghost_type)));

      Vector<UInt> conn(2);
      conn(0) = 0;

      bool found = false;
      GhostType ghost_type1 = _not_ghost;

      for (ghost_type_t::iterator git = ghost_type_t::begin();  git != ghost_type_t::end(); ++git) {
        GhostType ghost_type2 = *git;
        UInt existing_pairs_num = ghost_type2 - _not_ghost;
        pair_type::iterator first_pair_types = existing_pairs[existing_pairs_num].begin();
        pair_type::iterator last_pair_types = existing_pairs[existing_pairs_num].end();

        // Compute the weights
        for (; first_pair_types != last_pair_types; ++first_pair_types) {
          ElementType type1 = first_pair_types->first;
          ElementType type2 = first_pair_types->second;

          const Array<UInt> & elem_mat_1 = element_filter(type1, ghost_type1);
          const Array<UInt> & elem_mat_2 = element_filter(type2, ghost_type2);

          const Array<UInt> & pairs = pair_list(type1, ghost_type1)(type2, ghost_type2);
          const Array<Real> & weights = pair_weight(type1, ghost_type1)(type2, ghost_type2);

          UInt nb_quad1 = fem.getNbQuadraturePoints(type1, ghost_type1);
          UInt nb_quad2 = fem.getNbQuadraturePoints(type2, ghost_type2);

          Array<UInt>::const_iterator< Vector<UInt> > first_pair = pairs.begin(2);
          Array<UInt>::const_iterator< Vector<UInt> > last_pair  = pairs.end(2);
          Array<Real>::const_vector_iterator pair_w = weights.begin(2);

          for(;first_pair != last_pair; ++first_pair, ++pair_w) {
            UInt _q1 = (*first_pair)(0);
            UInt _q2 = (*first_pair)(1);
            QuadraturePoint q1(type1, elem_mat_1(_q1 / nb_quad1), _q1 % nb_quad1, ghost_type1);
            QuadraturePoint q2(type2, elem_mat_2(_q2 / nb_quad2), _q2 % nb_quad2, ghost_type2);

            if(el == (Element) q1 || el == (Element) q2) {
              QuadraturePoint q;
	      Real w1;
              UInt ele;
              //              Real w2;
              if(el == (Element) q1) {
		q = q2;
		w1 = (*pair_w)(0);
                ele = q.element * nb_quad2;
                //                w2 = (*pair_w)(1);
	      } else {
		q = q1;
		w1 = (*pair_w)(1);
                ele = q.element * nb_quad1;
                //                w2 = (*pair_w)(0);
	      }

              found = true;

              fem.getNbQuadraturePoints(q.type, q.ghost_type);
              quad_it = quadrature_points_coordinates(q.type, q.ghost_type).begin(spatial_dimension);
              const Vector<Real> & quad_coord = quad_it[q.element * nb_quad_per_elem];
              std::map<Vector<Real>, UInt>::iterator nit = numbering.find(quad_coord);
              if(nit != numbering.end()) {
                conn(1) = nit->second;
                if(conn(1) == 0) {
                  weight(0) = w1;
                }
              } else {
                conn(1) = counter++;
                numbering[quad_coord] = conn(1);
                nodes.push_back(quad_coord);
                g_notg.push_back(q.ghost_type == _ghost);
                dist.push_back(quad.distance(quad_coord));
                weight.push_back(w1);
                damage.push_back((getArray("damage", q.type, q.ghost_type)(q.global_num)));
                jac.push_back(fem.getIntegratorInterface().getJacobians(q.type, q.ghost_type)(ele));
              }
              connectivity.push_back(conn);
            }
          }
        }
      }
      if(found) {
        std::stringstream sstr;
        sstr << "neigh_mesh" << el;
        out << "./neighbors" << sstr.str() << " dumped!";
        DumperParaview dumper(sstr.str(), "./neighbors", false);
        dumper.registerMesh(mesh);
        dumper.registerField("ghost", new DumperIOHelper::NodalField<bool>(g_notg));
        dumper.registerField("distance", new DumperIOHelper::NodalField<Real>(dist));
        dumper.registerField("weight", new DumperIOHelper::NodalField<Real>(weight));
        dumper.registerField("damage", new DumperIOHelper::NodalField<Real>(damage));
        dumper.registerField("jacobian", new DumperIOHelper::NodalField<Real>(jac));
        dumper.dump();
      }

      return out.str();
    });
#endif
  // <<<<<< DEBUG CODE <<<<<< //


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
void MaterialNonLocal<spatial_dimension, WeightFunction>::cleanupExtraGhostElement(const ElementTypeMap<UInt> & nb_ghost_protected) {
  AKANTU_DEBUG_IN();

  // Create list of element to keep
  std::set<Element> relevant_ghost_element;

  pair_type::const_iterator first_pair_types = existing_pairs[1].begin();
  pair_type::const_iterator last_pair_types = existing_pairs[1].end();
  for (; first_pair_types != last_pair_types; ++first_pair_types) {
    ElementType type2 = first_pair_types->second;
    GhostType ghost_type2 = _ghost;
    UInt nb_quad2 = this->model->getFEEngine().getNbQuadraturePoints(type2);
    Array<UInt> & elem_filter = element_filter(type2, ghost_type2);

    const Array<UInt> & pairs =
      pair_list(first_pair_types->first, _not_ghost)(first_pair_types->second, ghost_type2);
    Array<UInt>::const_iterator< Vector<UInt> > first_pair = pairs.begin(2);
    Array<UInt>::const_iterator< Vector<UInt> > last_pair  = pairs.end(2);
    for(;first_pair != last_pair; ++first_pair) {
      UInt _q2 = (*first_pair)(1);
      QuadraturePoint q2(type2, elem_filter(_q2 / nb_quad2), _q2 % nb_quad2, ghost_type2);
      relevant_ghost_element.insert(q2);
    }
  }

  // Create list of element to remove and new numbering for element to keep
  Mesh & mesh = this->model->getFEEngine().getMesh();
  std::set<Element> ghost_to_erase;

  Mesh::type_iterator it = mesh.firstType(spatial_dimension, _ghost);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, _ghost);

  RemovedElementsEvent remove_elem(mesh);

  Element element;
  element.ghost_type = _ghost;
  for(; it != last_type; ++it) {
    element.type = *it;
    UInt nb_ghost_elem = mesh.getNbElement(*it, _ghost);
    UInt nb_ghost_elem_protected = 0;
    try {
      nb_ghost_elem_protected = nb_ghost_protected(*it, _ghost);
    } catch (...) {}

    if(!remove_elem.getNewNumbering().exists(*it, _ghost))
      remove_elem.getNewNumbering().alloc(nb_ghost_elem, 1, *it, _ghost);
    else remove_elem.getNewNumbering(*it, _ghost).resize(nb_ghost_elem);

    Array<UInt> & elem_filter = element_filter(*it, _ghost);
    Array<UInt> & new_numbering = remove_elem.getNewNumbering(*it, _ghost);
    UInt ng = 0;
    for (UInt g = 0; g < nb_ghost_elem; ++g) {
      element.element = elem_filter(g);
      if(element.element >= nb_ghost_elem_protected &&
         (std::find(relevant_ghost_element.begin(),
                    relevant_ghost_element.end(),
                    element) == relevant_ghost_element.end())) {
        ghost_to_erase.insert(element);
        remove_elem.getList().push_back(element);

        new_numbering(g) = UInt(-1);
      } else {
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

  computeQuadraturePointsCoordinates(quadrature_points_coordinates, _not_ghost);
  fillCellList(quadrature_points_coordinates, _not_ghost);

  is_creating_grid = true;
  SynchronizerRegistry & synch_registry = this->model->getSynchronizerRegistry();
  std::stringstream sstr; sstr << getID() << ":grid_synchronizer";
  grid_synchronizer = GridSynchronizer::createGridSynchronizer(mesh,
                                                               *spatial_grid,
                                                               sstr.str());
  synch_registry.registerSynchronizer(*grid_synchronizer, _gst_mnl_for_average);
  synch_registry.registerSynchronizer(*grid_synchronizer, _gst_mnl_weight);
  is_creating_grid = false;

#if not defined(AKANTU_NDEBUG)
  Mesh * mesh_tmp = NULL;
  if(AKANTU_DEBUG_TEST(dblDump)) {
    mesh_tmp = new Mesh(spatial_dimension, "mnl_grid");
    spatial_grid->saveAsMesh(*mesh_tmp);
    std::stringstream sstr_grid;
    StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
    Int prank = comm.whoAmI();

    sstr_grid << "material_non_local_grid_" << std::setfill('0') << std::setw(4) << prank << ".msh";
    mesh_tmp->write(sstr_grid.str());
    delete mesh_tmp;
  }
#endif

  this->computeQuadraturePointsCoordinates(quadrature_points_coordinates, _ghost);
  fillCellList(quadrature_points_coordinates, _ghost);

#if not defined(AKANTU_NDEBUG)
  if(AKANTU_DEBUG_TEST(dblDump)) {
    mesh_tmp = new Mesh(spatial_dimension, "mnl_grid");
    spatial_grid->saveAsMesh(*mesh_tmp);
    std::stringstream sstr_grid;
    StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
    Int prank = comm.whoAmI();

    sstr_grid.str(std::string());
    sstr_grid << "material_non_local_grid_ghost_" << std::setfill('0') << std::setw(4) << prank << ".msh";

    mesh_tmp->write(sstr_grid.str());
    delete mesh_tmp;
  }
#endif

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
void MaterialNonLocal<spatial_dimension, WeightFunction>::fillCellList(const ElementTypeMapArray<Real> & quadrature_points_coordinates,
                                                                       const GhostType & ghost_type) {
  Mesh & mesh = this->model->getFEEngine().getMesh();

  QuadraturePoint q;
  q.ghost_type = ghost_type;

  Mesh::type_iterator it        = mesh.firstType(spatial_dimension, ghost_type);
  Mesh::type_iterator last_type = mesh.lastType (spatial_dimension, ghost_type);
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

  Mesh & mesh = this->model->getFEEngine().getMesh();

  GhostType ghost_type = _not_ghost;
  QuadraturePoint quad_point;
  quad_point.ghost_type = ghost_type;

  // generate the pair of neighbor depending of the cell_list
  Mesh::type_iterator it        = mesh.firstType(spatial_dimension, ghost_type);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, ghost_type);
  for(; it != last_type; ++it) {
    // Preparing datas
    const Array<Real> & quads = quadrature_points_coordinates(*it, ghost_type);
    Array<Real>::const_vector_iterator first_quad = quads.begin(spatial_dimension);
    Array<Real>::const_vector_iterator last_quad  = quads.end(spatial_dimension);

    ElementTypeMapArray<UInt> & pairs = pair_list(ElementTypeMapArray<UInt>("pairs", getID(), memory_id),
                                          *it,
                                          ghost_type);

    ElementType current_element_type = _not_defined;
    GhostType current_ghost_type = _casper;
    UInt existing_pairs_num = 0;

    Array<UInt> * neighbors = NULL;
    Array<UInt>::const_iterator< Vector<UInt> > element_index_material_it;
    Array<Real>::const_vector_iterator quad_coord_it;

    UInt my_num_quad = 0;
    quad_point.type = *it;
    // loop over quad points
    for(;first_quad != last_quad; ++first_quad, ++my_num_quad) {
      SpatialGrid<QuadraturePoint>::CellID cell_id = spatial_grid->getCellID(*first_quad);

      SpatialGrid<QuadraturePoint>::neighbor_cells_iterator first_neigh_cell =
        spatial_grid->beginNeighborCells(cell_id);
      SpatialGrid<QuadraturePoint>::neighbor_cells_iterator last_neigh_cell =
        spatial_grid->endNeighborCells(cell_id);

      quad_point.element = my_num_quad / this->model->getFEEngine().getNbQuadraturePoints(*it,
                                                                                     quad_point.ghost_type);

      // loop over neighbors cells of the one containing the current quadrature
      // point
      for (; first_neigh_cell != last_neigh_cell; ++first_neigh_cell) {
        SpatialGrid<QuadraturePoint>::Cell::iterator first_neigh_quad =
          spatial_grid->beginCell(*first_neigh_cell);
        SpatialGrid<QuadraturePoint>::Cell::iterator last_neigh_quad =
          spatial_grid->endCell(*first_neigh_cell);

        // loop over the quadrature point in the current cell of the cell list
        for (;first_neigh_quad != last_neigh_quad; ++first_neigh_quad){
          QuadraturePoint quad = *first_neigh_quad;
          UInt nb_quad_per_elem =
            this->model->getFEEngine().getNbQuadraturePoints(quad.type,
                                                        quad.ghost_type);

          // little optimization to not search in the map at each quad points
          if(quad.type != current_element_type ||
             quad.ghost_type != current_ghost_type) {

            current_element_type = quad.type;
            current_ghost_type   = quad.ghost_type;
            existing_pairs_num = quad.ghost_type == _not_ghost ? 0 : 1;
            if(!pairs.exists(current_element_type, current_ghost_type)) {
              neighbors = &(pairs.alloc(0, 2,
                                        current_element_type,
                                        current_ghost_type));
            } else {
              neighbors = &(pairs(current_element_type,
                                  current_ghost_type));
            }
            existing_pairs[existing_pairs_num].insert(std::pair<ElementType,
                                                      ElementType>(*it,
                                                                   current_element_type));
            element_index_material_it = this->model->getElementIndexByMaterial(current_element_type,
                                                                               current_ghost_type).begin(2);
            quad_coord_it = quadrature_points_coordinates(current_element_type, current_ghost_type).begin(spatial_dimension);
          }

          const Vector<UInt> & el_mat = element_index_material_it[quad.element];
          UInt neigh_num_quad = el_mat(1) * nb_quad_per_elem + quad.num_point;
          const Vector<Real> & neigh_quad = quad_coord_it[neigh_num_quad];

          Real distance = first_quad->distance(neigh_quad);
          if(distance <= weight_func->getRadius() &&
             (quad.ghost_type == _ghost ||
              (quad.ghost_type == _not_ghost && my_num_quad <= neigh_num_quad))) { // storing only half lists
            UInt pair[2];
            pair[0] = my_num_quad;
            pair[1] = neigh_num_quad;
            neighbors->push_back(pair);

            // >>>>>> DEBUG CODE >>>>>> //
#if defined(AKANTU_DEBUG_TOOLS) && defined(AKANTU_CORE_CXX11)
            debug::element_manager.print
              (debug::_dm_material_non_local,
               [this, &first_quad, &neigh_quad, &distance,
                &quad_point, &quad](const Element & el)->std::string {
                std::stringstream out;
                if((Element) quad_point == el) {
                  out << " neigh1: " << quad << " -- " << *first_quad << " " << neigh_quad << " dist: " << distance;
                }
                if((Element) quad == el) {
                  out << " neigh2: " << quad_point << " -- " << *first_quad << " " << neigh_quad << " dist: " << distance;
                }
                return out.str();
              });
#endif
            // <<<<<< DEBUG CODE <<<<<< //
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

  GhostType ghost_type1;
  ghost_type1 = _not_ghost;

  InternalField<Real> quadrature_points_volumes("quadrature_points_volumes", *this);
  quadrature_points_volumes.initialize(1);

  const FEEngine & fem = this->model->getFEEngine();

  weight_func->updateInternals(quadrature_points_volumes);

  for(UInt gt = _not_ghost; gt <= _ghost; ++gt) {
    GhostType ghost_type2 = (GhostType) gt;
    UInt existing_pairs_num = gt - _not_ghost;
    pair_type::iterator first_pair_types = existing_pairs[existing_pairs_num].begin();
    pair_type::iterator last_pair_types = existing_pairs[existing_pairs_num].end();

    // Compute the weights
    for (; first_pair_types != last_pair_types; ++first_pair_types) {
      ElementType type1 = first_pair_types->first;
      ElementType type2 = first_pair_types->second;

      const Array<UInt> & pairs = pair_list(type1, ghost_type1)(type2, ghost_type2);

      std::string ghost_id = "";
      if (ghost_type1 == _ghost) ghost_id = ":ghost";

      ElementTypeMapArray<Real> & weights_type_1 = pair_weight(type1, ghost_type1);
      std::stringstream sstr; sstr << getID() << ":pair_weight:" << type1 << ghost_id;
      weights_type_1.setID(sstr.str());

      Array<Real> * tmp_weight = NULL;
      if(!weights_type_1.exists(type2, ghost_type2)) {
        tmp_weight = &(weights_type_1.alloc(0, 2, type2, ghost_type2));
      } else {
        tmp_weight = &(weights_type_1(type2, ghost_type2));
      }
      Array<Real> & weights = *tmp_weight;
      weights.resize(pairs.getSize());
      weights.clear();

      const Array<Real> & jacobians_1 = fem.getIntegratorInterface().getJacobians(type1, ghost_type1);
      const Array<Real> & jacobians_2 = fem.getIntegratorInterface().getJacobians(type2, ghost_type2);

      const Array<UInt> & elem_mat_1 = element_filter(type1, ghost_type1);
      const Array<UInt> & elem_mat_2 = element_filter(type2, ghost_type2);

      UInt nb_quad1 = fem.getNbQuadraturePoints(type1);
      UInt nb_quad2 = fem.getNbQuadraturePoints(type2);

      Array<Real> & quads_volumes1 = quadrature_points_volumes(type1, ghost_type1);
      Array<Real> & quads_volumes2 = quadrature_points_volumes(type2, ghost_type2);

      Array<Real>::const_vector_iterator iquads1;
      Array<Real>::const_vector_iterator iquads2;
      iquads1 = quadrature_points_coordinates(type1, ghost_type1).begin(spatial_dimension);
      iquads2 = quadrature_points_coordinates(type2, ghost_type2).begin(spatial_dimension);

      Array<UInt>::const_iterator< Vector<UInt> > first_pair = pairs.begin(2);
      Array<UInt>::const_iterator< Vector<UInt> > last_pair  = pairs.end(2);
      Array<Real>::vector_iterator weight  = weights.begin(2);

      this->weight_func->selectType(type1, ghost_type1, type2, ghost_type2);

      // Weight function
      for(;first_pair != last_pair; ++first_pair, ++weight) {
        UInt _q1 = (*first_pair)(0);
        UInt _q2 = (*first_pair)(1);
        const Vector<Real> & pos1 = iquads1[_q1];
        const Vector<Real> & pos2 = iquads2[_q2];
        QuadraturePoint q1(_q1 / nb_quad1, _q1 % nb_quad1, _q1, pos1, type1, ghost_type1);
        QuadraturePoint q2(_q2 / nb_quad2, _q2 % nb_quad2, _q2, pos2, type2, ghost_type2);

        Real r = pos1.distance(pos2);

        Real w2J2 = jacobians_2(elem_mat_2(q2.element)*nb_quad2 + q2.num_point);
        Real w = this->weight_func->operator()(r, q1, q2);
        (*weight)(0) = w2J2 * w;
        quads_volumes1(_q1) += (*weight)(0);

        if(ghost_type2 != _ghost && _q1 != _q2) {
          Real w1J1 = jacobians_1(elem_mat_1(q1.element)*nb_quad1 + q1.num_point);
          (*weight)(1) = w1J1 * this->weight_func->operator()(r, q2, q1);
          quads_volumes2(_q2) += (*weight)(1);
        } else
          (*weight)(1) = 0;
      }
    }
  }

  //normalize the weights
  for(UInt gt = _not_ghost; gt <= _ghost; ++gt) {
    GhostType ghost_type2 = (GhostType) gt;
    UInt existing_pairs_num = gt - _not_ghost;
    pair_type::iterator first_pair_types = existing_pairs[existing_pairs_num].begin();
    pair_type::iterator last_pair_types = existing_pairs[existing_pairs_num].end();

    first_pair_types = existing_pairs[existing_pairs_num].begin();
    for (; first_pair_types != last_pair_types; ++first_pair_types) {
      ElementType type1 = first_pair_types->first;
      ElementType type2 = first_pair_types->second;

      const Array<UInt> & pairs = pair_list(type1, ghost_type1)(type2, ghost_type2);
      Array<Real> & weights = pair_weight(type1, ghost_type1)(type2, ghost_type2);

      Array<Real> & quads_volumes1 = quadrature_points_volumes(type1, ghost_type1);
      Array<Real> & quads_volumes2 = quadrature_points_volumes(type2, ghost_type2);

      Array<UInt>::const_iterator< Vector<UInt> > first_pair = pairs.begin(2);
      Array<UInt>::const_iterator< Vector<UInt> > last_pair  = pairs.end(2);
      Array<Real>::vector_iterator weight  = weights.begin(2);

      for(;first_pair != last_pair; ++first_pair, ++weight) {
        UInt q1 = (*first_pair)(0);
        UInt q2 = (*first_pair)(1);

        (*weight)(0) *= 1. / quads_volumes1(q1);
        if(ghost_type2 != _ghost) (*weight)(1) *= 1. / quads_volumes2(q2);
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
template<typename T>
void MaterialNonLocal<spatial_dimension, WeightFunction>::weightedAvergageOnNeighbours(const ElementTypeMapArray<T> & to_accumulate,
                                                                                       ElementTypeMapArray<T> & accumulated,
                                                                                       UInt nb_degree_of_freedom,
                                                                                       GhostType ghost_type2) const {
  AKANTU_DEBUG_IN();

  UInt existing_pairs_num = 0;
  if (ghost_type2 == _ghost) existing_pairs_num = 1;

  pair_type::const_iterator first_pair_types = existing_pairs[existing_pairs_num].begin();
  pair_type::const_iterator last_pair_types = existing_pairs[existing_pairs_num].end();

  GhostType ghost_type1 = _not_ghost; // does not make sens the ghost vs ghost so this should always by _not_ghost

  for (; first_pair_types != last_pair_types; ++first_pair_types) {
    const Array<UInt> & pairs =
      pair_list(first_pair_types->first, ghost_type1)(first_pair_types->second, ghost_type2);
    const Array<Real> & weights =
      pair_weight(first_pair_types->first, ghost_type1)(first_pair_types->second, ghost_type2);

    const Array<T> & to_acc = to_accumulate(first_pair_types->second, ghost_type2);
    Array<T> & acc = accumulated(first_pair_types->first, ghost_type1);

    if(ghost_type2 == _not_ghost) acc.clear();

    Array<UInt>::const_iterator< Vector<UInt> > first_pair = pairs.begin(2);
    Array<UInt>::const_iterator< Vector<UInt> > last_pair  = pairs.end(2);
    Array<Real>::const_vector_iterator pair_w = weights.begin(2);

    for(;first_pair != last_pair; ++first_pair, ++pair_w) {
      UInt q1 = (*first_pair)(0);
      UInt q2 = (*first_pair)(1);
      for(UInt d = 0; d < nb_degree_of_freedom; ++d){
        acc(q1, d) += (*pair_w)(0) * to_acc(q2, d);
        if(ghost_type2 != _ghost) acc(q2, d) += (*pair_w)(1) * to_acc(q1, d);
      }
    }
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

      // >>>>>> DEBUG CODE >>>>>> //
#if defined(AKANTU_DEBUG_TOOLS)
#if defined(AKANTU_CORE_CXX11)
      debug::element_manager.print(debug::_dm_material,
                                   [non_local_variable, &element_filter](const Element & el)->std::string {
                                     std::stringstream out;
                                     if(element_filter.exists(el.type, el.ghost_type)) {
                                       UInt pos = element_filter(el.type, el.ghost_type).find(el.element);
                                       if(pos != UInt(-1)) {
                                         out << (*non_local_variable.local)(el.type, el.ghost_type).getID()
                                             << ": loc "
                                             << (*non_local_variable.local)(el.type, el.ghost_type).begin(non_local_variable.nb_component)[pos]
                                             << " - non loc "
                                             << (*non_local_variable.non_local)(el.type, el.ghost_type).begin(non_local_variable.nb_component)[pos];
                                       }
                                     }
                                     return out.str();
                                   });
#else
      debug::element_manager.printData(debug::_dm_material, "MaterialNonLocal: computeAllNonLocalStress",
                                       *non_local_variable.local, this->element_filter);
      debug::element_manager.printData(debug::_dm_material, "MaterialNonLocal: computeAllNonLocalStress",
                                       *non_local_variable.non_local, this->element_filter);
#endif
#endif
      // <<<<<< DEBUG CODE <<<<<< //
    }

    computeNonLocalStresses(_not_ghost);

  // >>>>>> DEBUG CODE >>>>>> //
#if defined(AKANTU_DEBUG_TOOLS) && defined(AKANTU_CORE_CXX11)
  ElementTypeMapArray<Real> quadrature_points_coordinates("quadrature_points_coordinates", id);
  Mesh & mesh = this->model->getFEEngine().getMesh();
  mesh.initElementTypeMapArray(quadrature_points_coordinates, spatial_dimension, spatial_dimension);
  computeQuadraturePointsCoordinates(quadrature_points_coordinates, _not_ghost);
  computeQuadraturePointsCoordinates(quadrature_points_coordinates, _ghost);

  debug::element_manager.print
    (debug::_dm_material_non_local,
     [this, &quadrature_points_coordinates](const Element & el)->std::string {
      static UInt step = 0;
      std::stringstream out;

      GhostType ghost_type1 = _not_ghost;
      FEEngine & fem = this->model->getFEEngine();
      Mesh & mesh = this->model->getMesh();

      std::ofstream quad_out;
      std::stringstream sstro;
      sstro << "neigh_vals_" << el.element << ".csv";
      if(step == 0) {
        quad_out.open(sstro.str());
        quad_out << "#id,step,gt,realid,w1,w2";
        for (UInt i = 0; i < spatial_dimension; ++i) {
          std::stringstream sstr; sstr << ",x" << i ;
          quad_out << sstr.str();
        }
        quad_out << ",dam,jac";

        typename std::map<ID, NonLocalVariable>::iterator it = non_local_variables.begin();
        typename std::map<ID, NonLocalVariable>::iterator end = non_local_variables.end();
        for(;it != end; ++it) {
          NonLocalVariable & non_local_variable = it->second;
          for (UInt i = 0; i < non_local_variable.nb_component; ++i) {
            std::stringstream sstr; sstr << "," << non_local_variable.local->getID() << i;
            quad_out << sstr.str();
          }
          for (UInt i = 0; i < non_local_variable.nb_component; ++i) {
            std::stringstream sstr; sstr << "," << non_local_variable.non_local->getID() << i;
            quad_out << sstr.str();
          }
        }
        quad_out << std::endl;
      } else {
        quad_out.open(sstro.str(), std::ios_base::app);
      }
      quad_out.precision(16);

      for (ghost_type_t::iterator git = ghost_type_t::begin();  git != ghost_type_t::end(); ++git) {
        GhostType ghost_type2 = *git;
        UInt existing_pairs_num = ghost_type2 - _not_ghost;
        pair_type::iterator first_pair_types = existing_pairs[existing_pairs_num].begin();
        pair_type::iterator last_pair_types = existing_pairs[existing_pairs_num].end();

        // Compute the weights
        for (; first_pair_types != last_pair_types; ++first_pair_types) {
          ElementType type1 = first_pair_types->first;
          ElementType type2 = first_pair_types->second;

          const Array<UInt> & elem_mat_1 = element_filter(type1, ghost_type1);
          const Array<UInt> & elem_mat_2 = element_filter(type2, ghost_type2);

          UInt nb_quad1 = fem.getNbQuadraturePoints(type1, ghost_type1);
          UInt nb_quad2 = fem.getNbQuadraturePoints(type2, ghost_type2);

          const Array<UInt> & pairs = pair_list(type1, ghost_type1)(type2, ghost_type2);
          const Array<Real> & weights = pair_weight(type1, ghost_type1)(type2, ghost_type2);

          Array<UInt>::const_iterator< Vector<UInt> > first_pair = pairs.begin(2);
          Array<UInt>::const_iterator< Vector<UInt> > last_pair  = pairs.end(2);
          Array<Real>::const_vector_iterator pair_w = weights.begin(2);

          for(;first_pair != last_pair; ++first_pair, ++pair_w) {
            UInt _q1 = (*first_pair)(0);
            UInt _q2 = (*first_pair)(1);

            QuadraturePoint q1(type1, elem_mat_1(_q1 / nb_quad1), _q1 % nb_quad1, ghost_type1);
            QuadraturePoint q2(type2, elem_mat_2(_q2 / nb_quad2), _q2 % nb_quad2, ghost_type2);
            q1.global_num = _q1;
            q2.global_num = _q2;

            if(el == (Element) q1 || el == (Element) q2) {
              QuadraturePoint q;
	      Real w1, w2;
              UInt nb_quad;
              if(el == (Element) q1) {
		q = q2;
                nb_quad = nb_quad2;
		w1 = (*pair_w)(0);
		w2 = (*pair_w)(1);
	      } else {
		q = q1;
                nb_quad = nb_quad1;
		w1 = (*pair_w)(1);
		w2 = (*pair_w)(0);
	      }

              quad_out << ((q.ghost_type == _ghost ? mesh.getNbElement(q.type, q.ghost_type) : 0) + q.global_num)
		       << "," << step << "," << q.ghost_type << "," << q.global_num;
              quad_out << "," << w1 << "," << w2;
              Array<Real>::const_vector_iterator iquads =
                quadrature_points_coordinates(q.type, q.ghost_type).begin(spatial_dimension);
              const Vector<Real> & coord = iquads[q.global_num];
              for(UInt i(0); i < coord.size(); ++i) {
                quad_out << "," << coord(i);
              }

              typename std::map<ID, NonLocalVariable>::iterator nlit = non_local_variables.begin();
              typename std::map<ID, NonLocalVariable>::iterator nlend = non_local_variables.end();
              for(;nlit != nlend; ++nlit) {
                NonLocalVariable & non_local_variable = nlit->second;
                const Array<Real> & local = (*non_local_variable.local)(q.type, q.ghost_type);
                for(UInt i(0); i < local.getNbComponent(); ++i) {
                  quad_out << "," << local(q.global_num, i);
                }
                const Array<Real> & non_local = (*non_local_variable.non_local)(q.type, q.ghost_type);
                for(UInt i(0); i < non_local.getNbComponent(); ++i) {
                  quad_out << "," << non_local(q.global_num, i);
                }
              }
              Real d = getArray("damage", q.type, q.ghost_type)(q.global_num);
              Real j = fem.getIntegratorInterface().getJacobians(q.type, q.ghost_type)(q.element*nb_quad + q.num_point);
              quad_out << "," << d << "," << j;

              quad_out << std::endl;
            }
          }
        }
      }
      step++;
      return out.str();
    });
#endif
  // <<<<<< DEBUG CODE <<<<<< //
  }

  // >>>>>> DEBUG CODE >>>>>> //
#if defined(AKANTU_DEBUG_TOOLS)
#if defined(AKANTU_CORE_CXX11)
  debug::element_manager.print(debug::_dm_material_damage,
                               [ghost_type, this](const Element & el)->std::string {
                                 std::stringstream out;
                                 if(el.ghost_type == ghost_type && element_filter.exists(el.type, ghost_type)) {
                                   UInt pos = element_filter(el.type, el.ghost_type).find(el.element);
                                   if(pos != UInt(-1)) {
                                     Real d = getArray("damage", el.type, el.ghost_type)(pos);
                                     out << " damage: " << d;
                                   }
                                 }
                                 return out.str();
                               });
#endif
#endif
  // <<<<<< DEBUG CODE <<<<<< //

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

  GhostType ghost_type1;
  ghost_type1 = _not_ghost;

  for(UInt gt = _not_ghost; gt <= _ghost; ++gt) {
    GhostType ghost_type2 = (GhostType) gt;
    UInt existing_pairs_num = gt - _not_ghost;

    pair_type::const_iterator first_pair_types = existing_pairs[existing_pairs_num].begin();
    pair_type::const_iterator last_pair_types = existing_pairs[existing_pairs_num].end();

    for (; first_pair_types != last_pair_types; ++first_pair_types) {
      const Array<UInt> & pairs =
        (*pair_list(first_pair_types->first, ghost_type1))(first_pair_types->second, ghost_type2);
      const Array<Real> & weights =
        (*pair_weight(first_pair_types->first, ghost_type1))(first_pair_types->second, ghost_type2);

      pout << "Types : " << first_pair_types->first << " (" << ghost_type1 << ") - " << first_pair_types->second << " (" << ghost_type2 << ")" << std::endl;

      Array<UInt>::const_iterator< Vector<UInt> > first_pair = pairs.begin(2);
      Array<UInt>::const_iterator< Vector<UInt> > last_pair  = pairs.end(2);
      Array<Real>::const_vector_iterator pair_w = weights.begin(2);

      for(;first_pair != last_pair; ++first_pair, ++pair_w) {
        UInt q1 = (*first_pair)(0);
        UInt q2 = (*first_pair)(1);
        pout << q1 << " " << q2 << " "<< (*pair_w)(0) << " " << (*pair_w)(1) << std::endl;
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
void MaterialNonLocal<spatial_dimension, WeightFunction>::neighbourhoodStatistics(const std::string & filename) const {
  std::ofstream pout;
  pout.open(filename.c_str());

  const Mesh & mesh = this->model->getFEEngine().getMesh();

  GhostType ghost_type1;
  ghost_type1 = _not_ghost;

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  Int prank = comm.whoAmI();

  InternalField<UInt> nb_neighbors("nb_neighbours", *const_cast<MaterialNonLocal *>(this));
  nb_neighbors.initialize(1);

  for(UInt gt = _not_ghost; gt <= _ghost; ++gt) {
    GhostType ghost_type2 = (GhostType) gt;
    UInt existing_pairs_num = gt - _not_ghost;

    pair_type::const_iterator first_pair_types = existing_pairs[existing_pairs_num].begin();
    pair_type::const_iterator last_pair_types = existing_pairs[existing_pairs_num].end();

    for (; first_pair_types != last_pair_types; ++first_pair_types) {
      const Array<UInt> & pairs =
        pair_list(first_pair_types->first, ghost_type1)(first_pair_types->second, ghost_type2);
      if(prank == 0) {
        pout << ghost_type2 << " ";
        pout << "Types : " << first_pair_types->first << " " << first_pair_types->second << std::endl;
      }
      Array<UInt>::const_iterator< Vector<UInt> > first_pair = pairs.begin(2);
      Array<UInt>::const_iterator< Vector<UInt> > last_pair  = pairs.end(2);
      Array<UInt> & nb_neigh_1 = nb_neighbors(first_pair_types->first, ghost_type1);
      Array<UInt> & nb_neigh_2 = nb_neighbors(first_pair_types->second, ghost_type2);
      for(;first_pair != last_pair; ++first_pair) {
        UInt q1 = (*first_pair)(0);
        UInt q2 = (*first_pair)(1);
        ++(nb_neigh_1(q1));
        if(q1 != q2) ++(nb_neigh_2(q2));
      }
    }

    Mesh::type_iterator it        = mesh.firstType(spatial_dimension, ghost_type1);
    Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, ghost_type1);
    UInt nb_quads = 0;
    Real sum_nb_neig = 0;
    UInt max_nb_neig = 0;
    UInt min_nb_neig = std::numeric_limits<UInt>::max();
    for(; it != last_type; ++it) {
      Array<UInt> & nb_neighor = nb_neighbors(*it, ghost_type1);
      Array<UInt>::iterator<UInt> nb_neigh = nb_neighor.begin();
      Array<UInt>::iterator<UInt> end_neigh  = nb_neighor.end();

      for (; nb_neigh != end_neigh; ++nb_neigh, ++nb_quads) {
        UInt nb = *nb_neigh;
        sum_nb_neig += nb;
        max_nb_neig = std::max(max_nb_neig, nb);
        min_nb_neig = std::min(min_nb_neig, nb);
      }
    }


    comm.allReduce(&nb_quads,    1, _so_sum);
    comm.allReduce(&sum_nb_neig, 1, _so_sum);
    comm.allReduce(&max_nb_neig, 1, _so_max);
    comm.allReduce(&min_nb_neig, 1, _so_min);

    if(prank == 0) {
      pout << ghost_type2 << " ";
      pout << "Nb quadrature points: " << nb_quads << std::endl;

      Real mean_nb_neig = sum_nb_neig / Real(nb_quads);
      pout << ghost_type2 << " ";
      pout << "Average nb neighbors: " << mean_nb_neig << "(" << sum_nb_neig << ")" << std::endl;

      pout << ghost_type2 << " ";
      pout << "Max nb neighbors:     " << max_nb_neig << std::endl;

      pout << ghost_type2 << " ";
      pout << "Min nb neighbors:     " << min_nb_neig << std::endl;
    }
  }
  pout.close();
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
// template<UInt spatial_dimension, template <UInt> class WeightFunction>
// inline void MaterialNonLocal<spatial_dimension, WeightFunction>::onElementsAdded(const Array<Element> & element_list) {
//   AKANTU_DEBUG_IN();

//   Material::onElementsAdded(element_list, event);

//   AKANTU_DEBUG_OUT();
// }

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeightFunction>
inline void MaterialNonLocal<spatial_dimension, WeightFunction>::onElementsRemoved(const Array<Element> & element_list,
                                                                                   const ElementTypeMapArray<UInt> & new_numbering,
                                                                                   __attribute__((unused)) const RemovedElementsEvent & event) {
  AKANTU_DEBUG_IN();

  Material::onElementsRemoved(element_list, new_numbering, event);

  pair_type::const_iterator first_pair_types = existing_pairs[1].begin();
  pair_type::const_iterator last_pair_types = existing_pairs[1].end();

  // Renumber element to keep
  for (; first_pair_types != last_pair_types; ++first_pair_types) {
    ElementType type2 = first_pair_types->second;
    GhostType ghost_type2 = _ghost;
    UInt nb_quad2 = this->model->getFEEngine().getNbQuadraturePoints(type2);

    Array<UInt> & pairs =
      pair_list(first_pair_types->first, _not_ghost)(first_pair_types->second, ghost_type2);
    Array<UInt>::iterator< Vector<UInt> > first_pair = pairs.begin(2);
    Array<UInt>::iterator< Vector<UInt> > last_pair  = pairs.end(2);
    for(;first_pair != last_pair; ++first_pair) {
      UInt _q2 = (*first_pair)(1);
      const Array<UInt> & renumbering = new_numbering(type2, ghost_type2);
      UInt el = _q2 / nb_quad2;
      UInt new_el = renumbering(el);
      AKANTU_DEBUG_ASSERT(new_el != UInt(-1), "A local quad as been removed instead f just renumbered");
      (*first_pair)(1) = new_el * nb_quad2 + _q2 % nb_quad2;
    }
  }

  AKANTU_DEBUG_OUT();
}
