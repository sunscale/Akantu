/**
 * @file   non_local_manager.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Apr 13 2012
 * @date last modification: Wed Dec 16 2015
 *
 * @brief  Implementation of non-local manager
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
#include "non_local_manager.hh"
#include "base_weight_function.hh"
#include "material_non_local.hh"
#include "non_local_neighborhood.hh"
/* -------------------------------------------------------------------------- */
#include <numeric>
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
NonLocalManager::NonLocalManager(SolidMechanicsModel & model, const ID & id,
                                 const MemoryID & memory_id)
    : Memory(id, memory_id), Parsable(_st_neighborhoods, id), model(model),
      quad_positions("quad_positions", id, memory_id),
      volumes("volumes", id, memory_id),
      spatial_dimension(this->model.getSpatialDimension()),
      compute_stress_calls(0), dummy_registry(NULL), dummy_grid(NULL) {
  Mesh & mesh = this->model.getMesh();
  mesh.registerEventHandler(*this);

  /// initialize the element type map array
  /// it will be resized to nb_quad * nb_element during the computation of
  /// coords
  mesh.initElementTypeMapArray(quad_positions, spatial_dimension,
                               spatial_dimension, false, _ek_regular, true);
  this->initElementTypeMap(1, volumes, this->model.getFEEngine());

  /// parse the neighborhood information from the input file
  const Parser & parser = getStaticParser();

  /// iterate over all the non-local sections and store them in a map
  std::pair<Parser::const_section_iterator, Parser::const_section_iterator>
      weight_sect = parser.getSubSections(_st_non_local);
  Parser::const_section_iterator it = weight_sect.first;
  for (; it != weight_sect.second; ++it) {
    const ParserSection & section = *it;
    ID name = section.getName();
    this->weight_function_types[name] = section;
  }

  this->dummy_registry = new SynchronizerRegistry(this->dummy_accessor);
}

/* -------------------------------------------------------------------------- */
NonLocalManager::~NonLocalManager() {

  /// delete neighborhoods
  NeighborhoodMap::iterator it;
  for (it = neighborhoods.begin(); it != neighborhoods.end(); ++it) {
    if (it->second)
      delete it->second;
  }

  /// delete non-local variables
  std::map<ID, NonLocalVariable *>::iterator it_variables;
  for (it_variables = non_local_variables.begin();
       it_variables != non_local_variables.end(); ++it_variables) {
    if (it_variables->second)
      delete it_variables->second;
  }

  std::map<ID, ElementTypeMapReal *>::iterator it_internals;
  for (it_internals = weight_function_internals.begin();
       it_internals != weight_function_internals.end(); ++it_internals) {
    if (it_internals->second)
      delete it_internals->second;
  }

  std::map<ID, GridSynchronizer *>::iterator grid_synch_it;
  for (grid_synch_it = dummy_synchronizers.begin();
       grid_synch_it != dummy_synchronizers.end(); ++grid_synch_it) {
    if (grid_synch_it->second)
      delete grid_synch_it->second;
  }
  /// delete all objects related to the dummy synchronizers
  delete dummy_registry;
  delete dummy_grid;
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::setJacobians(const FEEngine & fe_engine,
                                   const ElementKind & kind) {
  Mesh & mesh = this->model.getMesh();
  for (UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType)g;
    Mesh::type_iterator it = mesh.firstType(spatial_dimension, gt, kind);
    Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, gt, kind);
    for (; it != last_type; ++it) {
      jacobians(*it, gt) =
          &fe_engine.getIntegratorInterface().getJacobians(*it, gt);
    }
  }
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::createNeighborhood(const ID & weight_func,
                                         const ID & neighborhood_id) {

  AKANTU_DEBUG_IN();

  const ParserSection & section = this->weight_function_types[weight_func];
  const ID weight_func_type = section.getOption();
  /// create new neighborhood for given ID
  std::stringstream sstr;
  sstr << id << ":neighborhood:" << neighborhood_id;

  if (weight_func_type == "base_wf")
    neighborhoods[neighborhood_id] =
        new NonLocalNeighborhood<BaseWeightFunction>(
            *this, this->quad_positions, sstr.str());
  else if (weight_func_type == "remove_wf")
    neighborhoods[neighborhood_id] =
        new NonLocalNeighborhood<RemoveDamagedWeightFunction>(
            *this, this->quad_positions, sstr.str());
  else if (weight_func_type == "stress_wf")
    neighborhoods[neighborhood_id] =
        new NonLocalNeighborhood<StressBasedWeightFunction>(
            *this, this->quad_positions, sstr.str());
  else if (weight_func_type == "damage_wf")
    neighborhoods[neighborhood_id] =
        new NonLocalNeighborhood<DamagedWeightFunction>(
            *this, this->quad_positions, sstr.str());
  else
    AKANTU_EXCEPTION("error in weight function type provided in material file");

  neighborhoods[neighborhood_id]->parseSection(section);
  neighborhoods[neighborhood_id]->initNeighborhood();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::createNeighborhoodSynchronizers() {
  /// exchange all the neighborhood IDs, so that every proc knows how many
  /// neighborhoods exist globally
  /// First: Compute locally the maximum ID size
  UInt max_id_size = 0;
  UInt current_size = 0;
  NeighborhoodMap::const_iterator it;
  for (it = neighborhoods.begin(); it != neighborhoods.end(); ++it) {
    current_size = it->first.size();
    if (current_size > max_id_size)
      max_id_size = current_size;
  }

  /// get the global maximum ID size on each proc
  StaticCommunicator & static_communicator =
      akantu::StaticCommunicator::getStaticCommunicator();
  static_communicator.allReduce(&max_id_size, 1, _so_max);

  /// get the rank for this proc and the total nb proc
  UInt prank = static_communicator.whoAmI();
  UInt psize = static_communicator.getNbProc();

  /// exchange the number of neighborhoods on each proc
  Array<Int> nb_neighborhoods_per_proc(psize);
  nb_neighborhoods_per_proc(prank) = neighborhoods.size();
  static_communicator.allGather(nb_neighborhoods_per_proc.storage(), 1);

  /// compute the total number of neighborhoods
  UInt nb_neighborhoods_global = std::accumulate(
      nb_neighborhoods_per_proc.begin(), nb_neighborhoods_per_proc.end(), 0);

  /// allocate an array of chars to store the names of all neighborhoods
  Array<char> buffer(nb_neighborhoods_global, max_id_size);

  /// starting index on this proc
  UInt starting_index =
      std::accumulate(nb_neighborhoods_per_proc.begin(),
                      nb_neighborhoods_per_proc.begin() + prank, 0);

  it = neighborhoods.begin();
  /// store the names of local neighborhoods in the buffer
  for (UInt i = 0; i < neighborhoods.size(); ++i, ++it) {
    UInt c = 0;
    for (; c < it->first.size(); ++c)
      buffer(i + starting_index, c) = it->first[c];

    for (; c < max_id_size; ++c)
      buffer(i + starting_index, c) = char(0);
  }

  /// store the nb of data to send in the all gather
  Array<Int> buffer_size(nb_neighborhoods_per_proc);
  buffer_size *= max_id_size;
  /// exchange the names of all the neighborhoods with all procs
  static_communicator.allGatherV(buffer.storage(), buffer_size.storage());

  for (UInt i = 0; i < nb_neighborhoods_global; ++i) {
    std::stringstream neighborhood_id;
    for (UInt c = 0; c < max_id_size; ++c) {
      if (buffer(i, c) == char(0))
        break;
      neighborhood_id << buffer(i, c);
    }
    global_neighborhoods.insert(neighborhood_id.str());
  }

  /// this proc does not know all the neighborhoods -> create dummy
  /// grid so that this proc can participate in the all gather for
  /// detecting the overlap of neighborhoods this proc doesn't know
  Vector<Real> grid_center(this->spatial_dimension);
  for (UInt s = 0; s < this->spatial_dimension; ++s)
    grid_center(s) = std::numeric_limits<Real>::max();
  dummy_grid =
      new SpatialGrid<IntegrationPoint>(spatial_dimension, 0., grid_center);
  std::set<SynchronizationTag> tags;
  tags.insert(_gst_mnl_for_average);
  tags.insert(_gst_mnl_weight);

  std::set<ID>::const_iterator global_neighborhoods_it =
      global_neighborhoods.begin();
  for (; global_neighborhoods_it != global_neighborhoods.end();
       ++global_neighborhoods_it) {
    it = neighborhoods.find(*global_neighborhoods_it);
    if (it != neighborhoods.end()) {
      it->second->createGridSynchronizer();
    } else {
      ID neighborhood_name = *global_neighborhoods_it;
      std::stringstream sstr;
      sstr << getID() << ":" << neighborhood_name << ":grid_synchronizer";
      dummy_synchronizers[neighborhood_name] =
          GridSynchronizer::createGridSynchronizer(
              this->model.getMesh(), *dummy_grid, sstr.str(), dummy_registry,
              tags, 0, false);
    }
  }
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::flattenInternal(ElementTypeMapReal & internal_flat,
                                      const GhostType & ghost_type,
                                      const ElementKind & kind) {

  const ID field_name = internal_flat.getName();
  for (UInt m = 0; m < this->non_local_materials.size(); ++m) {
    Material & material = *(this->non_local_materials[m]);
    if (material.isInternal<Real>(field_name, kind))
      material.flattenInternal(field_name, internal_flat, ghost_type, kind);
  }
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::averageInternals(const GhostType & ghost_type) {
  /// update the weights of the weight function
  if (ghost_type == _not_ghost)
    this->computeWeights();

  /// loop over all neighborhoods and compute the non-local variables
  NeighborhoodMap::iterator neighborhood_it = neighborhoods.begin();
  NeighborhoodMap::iterator neighborhood_end = neighborhoods.end();
  for (; neighborhood_it != neighborhood_end; ++neighborhood_it) {
    /// loop over all the non-local variables of the given neighborhood
    std::map<ID, NonLocalVariable *>::iterator non_local_variable_it =
        non_local_variables.begin();
    std::map<ID, NonLocalVariable *>::iterator non_local_variable_end =
        non_local_variables.end();
    for (; non_local_variable_it != non_local_variable_end;
         ++non_local_variable_it) {
      NonLocalVariable * non_local_var = non_local_variable_it->second;
      neighborhood_it->second->weightedAverageOnNeighbours(
          non_local_var->local, non_local_var->non_local,
          non_local_var->nb_component, ghost_type);
    }
  }
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::init() {

  /// store the number of current ghost elements for each type in the mesh
  ElementTypeMap<UInt> nb_ghost_protected;
  Mesh & mesh = this->model.getMesh();
  Mesh::type_iterator it = mesh.firstType(spatial_dimension, _ghost);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, _ghost);
  for (; it != last_type; ++it)
    nb_ghost_protected(mesh.getNbElement(*it, _ghost), *it, _ghost);

  /// exchange the missing ghosts for the non-local neighborhoods
  this->createNeighborhoodSynchronizers();

  /// insert the ghost quadrature points of the non-local materials into the
  /// non-local neighborhoods
  for (UInt m = 0; m < this->non_local_materials.size(); ++m) {
    switch (spatial_dimension) {
    case 1:
      dynamic_cast<MaterialNonLocal<1> &>(*(this->non_local_materials[m]))
          .insertQuadsInNeighborhoods(_ghost);
      break;
    case 2:
      dynamic_cast<MaterialNonLocal<2> &>(*(this->non_local_materials[m]))
          .insertQuadsInNeighborhoods(_ghost);
      break;
    case 3:
      dynamic_cast<MaterialNonLocal<3> &>(*(this->non_local_materials[m]))
          .insertQuadsInNeighborhoods(_ghost);
      break;
    }
  }

  FEEngine & fee = this->model.getFEEngine();
  this->updatePairLists();
  /// cleanup the unneccessary ghost elements
  this->cleanupExtraGhostElements(nb_ghost_protected);
  this->setJacobians(fee, _ek_regular);
  this->initNonLocalVariables();
  this->computeWeights();
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::initNonLocalVariables() {
  /// loop over all the non-local variables
  std::map<ID, NonLocalVariable *>::iterator non_local_variable_it =
      non_local_variables.begin();
  std::map<ID, NonLocalVariable *>::iterator non_local_variable_end =
      non_local_variables.end();

  for (; non_local_variable_it != non_local_variable_end;
       ++non_local_variable_it) {
    NonLocalVariable & variable = *(non_local_variable_it->second);
    this->initElementTypeMap(variable.nb_component, variable.non_local,
                             this->model.getFEEngine());
  }
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::initElementTypeMap(UInt nb_component,
                                         ElementTypeMapReal & element_map,
                                         const FEEngine & fee,
                                         const ElementKind el_kind) {
  Mesh & mesh = this->model.getMesh();
  /// need to resize the arrays
  for (UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType)g;
    Mesh::type_iterator it = mesh.firstType(spatial_dimension, gt, el_kind);
    Mesh::type_iterator end = mesh.lastType(spatial_dimension, gt, el_kind);
    for (; it != end; ++it) {
      ElementType el_type = *it;
      UInt nb_element = mesh.getNbElement(*it, gt);
      UInt nb_quads = fee.getNbIntegrationPoints(*it, gt);
      if (!element_map.exists(el_type, gt)) {
        element_map.alloc(nb_element * nb_quads, nb_component, el_type, gt);
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::distributeInternals(ElementKind kind) {

  /// loop over all the non-local variables and copy back their values into the
  /// materials
  std::map<ID, NonLocalVariable *>::iterator non_local_variable_it =
      non_local_variables.begin();
  std::map<ID, NonLocalVariable *>::iterator non_local_variable_end =
      non_local_variables.end();
  for (; non_local_variable_it != non_local_variable_end;
       ++non_local_variable_it) {
    NonLocalVariable * non_local_var = non_local_variable_it->second;
    const ID field_name = non_local_var->non_local.getName();
    /// loop over all the materials
    for (UInt m = 0; m < this->non_local_materials.size(); ++m) {
      if (this->non_local_materials[m]->isInternal<Real>(field_name, kind))

        switch (spatial_dimension) {
        case 1:
          dynamic_cast<MaterialNonLocal<1> &>(*(this->non_local_materials[m]))
              .updateNonLocalInternals(non_local_var->non_local, field_name,
                                       non_local_var->nb_component);
          break;
        case 2:
          dynamic_cast<MaterialNonLocal<2> &>(*(this->non_local_materials[m]))
              .updateNonLocalInternals(non_local_var->non_local, field_name,
                                       non_local_var->nb_component);
          break;
        case 3:
          dynamic_cast<MaterialNonLocal<3> &>(*(this->non_local_materials[m]))
              .updateNonLocalInternals(non_local_var->non_local, field_name,
                                       non_local_var->nb_component);
          break;
        }
    }
  }
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::computeAllNonLocalStresses() {

  /// update the flattened version of the internals
  std::map<ID, NonLocalVariable *>::iterator non_local_variable_it =
      non_local_variables.begin();
  std::map<ID, NonLocalVariable *>::iterator non_local_variable_end =
      non_local_variables.end();

  for (; non_local_variable_it != non_local_variable_end;
       ++non_local_variable_it) {
    non_local_variable_it->second->local.clear();
    non_local_variable_it->second->non_local.clear();
    for (UInt gt = _not_ghost; gt <= _ghost; ++gt) {
      GhostType ghost_type = (GhostType)gt;
      this->flattenInternal(non_local_variable_it->second->local, ghost_type,
                            _ek_regular);
    }
  }

  this->volumes.clear();
  /// loop over all the neighborhoods and compute intiate the
  /// exchange of the non-local_variables
  // std::set<ID>::const_iterator global_neighborhood_it =
  // global_neighborhoods.begin();
  // NeighborhoodMap::iterator it;
  // for(; global_neighborhood_it != global_neighborhoods.end();
  // ++global_neighborhood_it) {
  //   it = neighborhoods.find(*global_neighborhood_it);
  //   if (it != neighborhoods.end())
  //     it->second->getSynchronizerRegistry().asynchronousSynchronize(_gst_mnl_for_average);
  //     else
  // 	dummy_synchronizers[*global_neighborhood_it]->asynchronousSynchronize(dummy_accessor,
  // _gst_mnl_for_average);
  // }

  NeighborhoodMap::iterator neighborhood_it = neighborhoods.begin();
  // NeighborhoodMap::iterator neighborhood_end = neighborhoods.end();

  for (; neighborhood_it != neighborhoods.end(); ++neighborhood_it) {
    neighborhood_it->second->getSynchronizerRegistry().asynchronousSynchronize(
        _gst_mnl_for_average);
  }

  this->averageInternals(_not_ghost);

  AKANTU_DEBUG_INFO("Wait distant non local stresses");

  /// loop over all the neighborhoods and block until all non-local
  /// variables have been exchanged
  // global_neighborhood_it = global_neighborhoods.begin();
  // for(; global_neighborhood_it != global_neighborhoods.end();
  // ++global_neighborhood_it) {
  //   it = neighborhoods.find(*global_neighborhood_it);
  //   if (it != neighborhoods.end())
  //     it->second->getSynchronizerRegistry().waitEndSynchronize(_gst_mnl_for_average);
  //   else
  //     dummy_synchronizers[*global_neighborhood_it]->waitEndSynchronize(dummy_accessor,
  //     _gst_mnl_for_average);
  // }

  neighborhood_it = neighborhoods.begin();
  for (; neighborhood_it != neighborhoods.end(); ++neighborhood_it) {
    neighborhood_it->second->getSynchronizerRegistry().waitEndSynchronize(
        _gst_mnl_for_average);
  }

  this->averageInternals(_ghost);

  /// copy the results in the materials
  this->distributeInternals(_ek_regular);
  /// loop over all the materials and update the weights
  for (UInt m = 0; m < this->non_local_materials.size(); ++m) {
    switch (spatial_dimension) {
    case 1:
      dynamic_cast<MaterialNonLocal<1> &>(*(this->non_local_materials[m]))
          .computeNonLocalStresses(_not_ghost);
      break;
    case 2:
      dynamic_cast<MaterialNonLocal<2> &>(*(this->non_local_materials[m]))
          .computeNonLocalStresses(_not_ghost);
      break;
    case 3:
      dynamic_cast<MaterialNonLocal<3> &>(*(this->non_local_materials[m]))
          .computeNonLocalStresses(_not_ghost);
      break;
    }
  }

  ++this->compute_stress_calls;
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::cleanupExtraGhostElements(
    ElementTypeMap<UInt> & nb_ghost_protected) {

  typedef std::set<Element> ElementSet;
  ElementSet relevant_ghost_elements;
  ElementSet to_keep_per_neighborhood;
  /// loop over all the neighborhoods and get their protected ghosts
  NeighborhoodMap::iterator neighborhood_it = neighborhoods.begin();
  NeighborhoodMap::iterator neighborhood_end = neighborhoods.end();
  for (; neighborhood_it != neighborhood_end; ++neighborhood_it) {
    neighborhood_it->second->cleanupExtraGhostElements(
        to_keep_per_neighborhood);
    ElementSet::const_iterator it = to_keep_per_neighborhood.begin();
    for (; it != to_keep_per_neighborhood.end(); ++it)
      relevant_ghost_elements.insert(*it);
    to_keep_per_neighborhood.clear();
  }

  /// remove all unneccessary ghosts from the mesh
  /// Create list of element to remove and new numbering for element to keep
  Mesh & mesh = this->model.getMesh();
  ElementSet ghost_to_erase;

  Mesh::type_iterator it = mesh.firstType(spatial_dimension, _ghost);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, _ghost);

  RemovedElementsEvent remove_elem(mesh);
  Element element;
  element.ghost_type = _ghost;

  for (; it != last_type; ++it) {
    element.type = *it;
    UInt nb_ghost_elem = mesh.getNbElement(*it, _ghost);
    UInt nb_ghost_elem_protected = 0;
    try {
      nb_ghost_elem_protected = nb_ghost_protected(*it, _ghost);
    } catch (...) {
    }

    if (!remove_elem.getNewNumbering().exists(*it, _ghost))
      remove_elem.getNewNumbering().alloc(nb_ghost_elem, 1, *it, _ghost);
    else
      remove_elem.getNewNumbering(*it, _ghost).resize(nb_ghost_elem);
    Array<UInt> & new_numbering = remove_elem.getNewNumbering(*it, _ghost);
    for (UInt g = 0; g < nb_ghost_elem; ++g) {
      element.element = g;
      if (element.element >= nb_ghost_elem_protected &&
          relevant_ghost_elements.find(element) ==
              relevant_ghost_elements.end()) {
        remove_elem.getList().push_back(element);
        new_numbering(element.element) = UInt(-1);
      }
    }
    /// renumber remaining ghosts
    UInt ng = 0;
    for (UInt g = 0; g < nb_ghost_elem; ++g) {
      if (new_numbering(g) != UInt(-1)) {
        new_numbering(g) = ng;
        ++ng;
      }
    }
  }

  it = mesh.firstType(spatial_dimension, _not_ghost);
  last_type = mesh.lastType(spatial_dimension, _not_ghost);
  for (; it != last_type; ++it) {
    UInt nb_elem = mesh.getNbElement(*it, _not_ghost);
    if (!remove_elem.getNewNumbering().exists(*it, _not_ghost))
      remove_elem.getNewNumbering().alloc(nb_elem, 1, *it, _not_ghost);
    Array<UInt> & new_numbering = remove_elem.getNewNumbering(*it, _not_ghost);
    for (UInt e = 0; e < nb_elem; ++e) {
      new_numbering(e) = e;
    }
  }
  mesh.sendEvent(remove_elem);
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::onElementsRemoved(
    const Array<Element> & element_list,
    const ElementTypeMapArray<UInt> & new_numbering,
    __attribute__((unused)) const RemovedElementsEvent & event) {

  FEEngine & fee = this->model.getFEEngine();
  this->removeIntegrationPointsFromMap(event.getNewNumbering(),
                                       spatial_dimension, quad_positions, fee,
                                       _ek_regular);
  this->removeIntegrationPointsFromMap(event.getNewNumbering(), 1, volumes, fee,
                                       _ek_regular);

  /// loop over all the neighborhoods and call onElementsRemoved
  std::set<ID>::const_iterator global_neighborhood_it =
      global_neighborhoods.begin();
  NeighborhoodMap::iterator it;
  for (; global_neighborhood_it != global_neighborhoods.end();
       ++global_neighborhood_it) {
    it = neighborhoods.find(*global_neighborhood_it);
    if (it != neighborhoods.end())
      it->second->onElementsRemoved(element_list, new_numbering, event);
    else
      dummy_synchronizers[*global_neighborhood_it]->onElementsRemoved(
          element_list, new_numbering, event);
  }
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::onElementsAdded(__attribute__((unused))
                                      const Array<Element> & element_list,
                                      __attribute__((unused))
                                      const NewElementsEvent & event) {
  this->resizeElementTypeMap(1, volumes, model.getFEEngine());
  this->resizeElementTypeMap(spatial_dimension, quad_positions,
                             model.getFEEngine());
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::resizeElementTypeMap(UInt nb_component,
                                           ElementTypeMapReal & element_map,
                                           const FEEngine & fee,
                                           const ElementKind el_kind) {
  Mesh & mesh = this->model.getMesh();
  for (UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType)g;
    Mesh::type_iterator it = mesh.firstType(spatial_dimension, gt, el_kind);
    Mesh::type_iterator end = mesh.lastType(spatial_dimension, gt, el_kind);
    for (; it != end; ++it) {
      UInt nb_element = mesh.getNbElement(*it, gt);
      UInt nb_quads = fee.getNbIntegrationPoints(*it, gt);
      if (!element_map.exists(*it, gt))
        element_map.alloc(nb_element * nb_quads, nb_component, *it, gt);
      else
        element_map(*it, gt).resize(nb_element * nb_quads);
    }
  }
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::removeIntegrationPointsFromMap(
    const ElementTypeMapArray<UInt> & new_numbering, UInt nb_component,
    ElementTypeMapReal & element_map, const FEEngine & fee,
    const ElementKind el_kind) {

  for (UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType)g;
    ElementTypeMapArray<UInt>::type_iterator it =
        new_numbering.firstType(_all_dimensions, gt, el_kind);
    ElementTypeMapArray<UInt>::type_iterator end =
        new_numbering.lastType(_all_dimensions, gt, el_kind);
    for (; it != end; ++it) {
      ElementType type = *it;
      if (element_map.exists(type, gt)) {
        const Array<UInt> & renumbering = new_numbering(type, gt);

        Array<Real> & vect = element_map(type, gt);

        UInt nb_quad_per_elem = fee.getNbIntegrationPoints(type, gt);

        Array<Real> tmp(renumbering.getSize() * nb_quad_per_elem, nb_component);

        AKANTU_DEBUG_ASSERT(
            tmp.getSize() == vect.getSize(),
            "Something strange append some mater was created from nowhere!!");

        AKANTU_DEBUG_ASSERT(
            tmp.getSize() == vect.getSize(),
            "Something strange append some mater was created or disappeared in "
                << vect.getID() << "(" << vect.getSize()
                << "!=" << tmp.getSize() << ") "
                                            "!!");

        UInt new_size = 0;
        for (UInt i = 0; i < renumbering.getSize(); ++i) {
          UInt new_i = renumbering(i);
          if (new_i != UInt(-1)) {
            memcpy(tmp.storage() + new_i * nb_component * nb_quad_per_elem,
                   vect.storage() + i * nb_component * nb_quad_per_elem,
                   nb_component * nb_quad_per_elem * sizeof(Real));
            ++new_size;
          }
        }
        tmp.resize(new_size * nb_quad_per_elem);
        vect.copy(tmp);
      }
    }
  }
}

__END_AKANTU__
