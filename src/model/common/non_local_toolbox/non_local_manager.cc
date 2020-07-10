/**
 * @file   non_local_manager.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Apr 13 2012
 * @date last modification: Tue Jan 16 2018
 *
 * @brief  Implementation of non-local manager
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
#include "non_local_manager.hh"
#include "grid_synchronizer.hh"
#include "model.hh"
#include "non_local_neighborhood.hh"
/* -------------------------------------------------------------------------- */
#include <numeric>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
NonLocalManager::NonLocalManager(Model & model,
                                 NonLocalManagerCallback & callback,
                                 const ID & id, const MemoryID & memory_id)
    : Memory(id, memory_id), Parsable(ParserType::_neighborhoods, id),
      spatial_dimension(model.getMesh().getSpatialDimension()), model(model),
      integration_points_positions("integration_points_positions", id,
                                   memory_id),
      volumes("volumes", id, memory_id), compute_stress_calls(0),
      dummy_registry(nullptr), dummy_grid(nullptr) {
  /// parse the neighborhood information from the input file
  const Parser & parser = getStaticParser();

  /// iterate over all the non-local sections and store them in a map
  std::pair<Parser::const_section_iterator, Parser::const_section_iterator>
      weight_sect = parser.getSubSections(ParserType::_non_local);
  Parser::const_section_iterator it = weight_sect.first;
  for (; it != weight_sect.second; ++it) {
    const ParserSection & section = *it;
    ID name = section.getName();
    this->weight_function_types[name] = section;
  }

  this->callback = &callback;
}

/* -------------------------------------------------------------------------- */
NonLocalManager::~NonLocalManager() = default;

/* -------------------------------------------------------------------------- */
void NonLocalManager::initialize() {
  volumes.initialize(this->model.getFEEngine(),
                     _spatial_dimension = spatial_dimension);

  AKANTU_DEBUG_ASSERT(this->callback,
                      "A callback should be registered prior to this call");
  this->callback->insertIntegrationPointsInNeighborhoods(_not_ghost);

  auto & mesh = this->model.getMesh();
  mesh.registerEventHandler(*this, _ehp_non_local_manager);

  /// store the number of current ghost elements for each type in the mesh
  // ElementTypeMap<UInt> nb_ghost_protected;
  // for (auto type : mesh.elementTypes(spatial_dimension, _ghost))
  //   nb_ghost_protected(mesh.getNbElement(type, _ghost), type, _ghost);

  /// exchange the missing ghosts for the non-local neighborhoods
  this->createNeighborhoodSynchronizers();

  /// insert the ghost quadrature points of the non-local materials into the
  /// non-local neighborhoods
  this->callback->insertIntegrationPointsInNeighborhoods(_ghost);

  FEEngine & fee = this->model.getFEEngine();
  this->updatePairLists();

  /// cleanup the unneccessary ghost elements
  this->cleanupExtraGhostElements(); // nb_ghost_protected);

  this->callback->initializeNonLocal();

  this->setJacobians(fee, _ek_regular);

  this->initNonLocalVariables();
  this->computeWeights();
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::setJacobians(const FEEngine & fe_engine,
                                   const ElementKind & kind) {
  Mesh & mesh = this->model.getMesh();
  for (auto ghost_type : ghost_types) {
    for (auto type : mesh.elementTypes(spatial_dimension, ghost_type, kind)) {
      jacobians(type, ghost_type) =
          &fe_engine.getIntegratorInterface().getJacobians(type, ghost_type);
    }
  }
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::createNeighborhood(const ID & weight_func,
                                         const ID & neighborhood_id) {

  AKANTU_DEBUG_IN();

  auto weight_func_it = this->weight_function_types.find(weight_func);
  AKANTU_DEBUG_ASSERT(weight_func_it != weight_function_types.end(),
                      "No info found in the input file for the weight_function "
                          << weight_func << " in the neighborhood "
                          << neighborhood_id);

  const ParserSection & section = weight_func_it->second;
  const ID weight_func_type = section.getOption();
  /// create new neighborhood for given ID
  std::stringstream sstr;
  sstr << id << ":neighborhood:" << neighborhood_id;

  if (weight_func_type == "base_wf")
    neighborhoods[neighborhood_id] =
        std::make_unique<NonLocalNeighborhood<BaseWeightFunction>>(
            *this, this->integration_points_positions, sstr.str());
#if defined(AKANTU_DAMAGE_NON_LOCAL)
  else if (weight_func_type == "remove_wf")
    neighborhoods[neighborhood_id] =
        std::make_unique<NonLocalNeighborhood<RemoveDamagedWeightFunction>>(
            *this, this->integration_points_positions, sstr.str());
  else if (weight_func_type == "stress_wf")
    neighborhoods[neighborhood_id] =
        std::make_unique<NonLocalNeighborhood<StressBasedWeightFunction>>(
            *this, this->integration_points_positions, sstr.str());
  else if (weight_func_type == "damage_wf")
    neighborhoods[neighborhood_id] =
        std::make_unique<NonLocalNeighborhood<DamagedWeightFunction>>(
            *this, this->integration_points_positions, sstr.str());
#endif
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
  const Communicator & static_communicator = model.getMesh().getCommunicator();
  static_communicator.allReduce(max_id_size, SynchronizerOperation::_max);

  /// get the rank for this proc and the total nb proc
  UInt prank = static_communicator.whoAmI();
  UInt psize = static_communicator.getNbProc();

  /// exchange the number of neighborhoods on each proc
  Array<Int> nb_neighborhoods_per_proc(psize);
  nb_neighborhoods_per_proc(prank) = neighborhoods.size();
  static_communicator.allGather(nb_neighborhoods_per_proc);

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
  static_communicator.allGatherV(buffer, buffer_size);

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
  Vector<Real> grid_center(this->spatial_dimension,
                           std::numeric_limits<Real>::max());
  Vector<Real> spacing(this->spatial_dimension, 0.);

  dummy_grid = std::make_unique<SpatialGrid<IntegrationPoint>>(
      this->spatial_dimension, spacing, grid_center);

  for (auto & neighborhood_id : global_neighborhoods) {
    it = neighborhoods.find(neighborhood_id);
    if (it != neighborhoods.end()) {
      it->second->createGridSynchronizer();
    } else {
      dummy_synchronizers[neighborhood_id] = std::make_unique<GridSynchronizer>(
          this->model.getMesh(), *dummy_grid,
          std::string(this->id + ":" + neighborhood_id + ":grid_synchronizer"),
          this->memory_id, false);
    }
  }
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::synchronize(DataAccessor<Element> & data_accessor,
                                  const SynchronizationTag & tag) {
  for (auto & neighborhood_id : global_neighborhoods) {
    auto it = neighborhoods.find(neighborhood_id);
    if (it != neighborhoods.end()) {
      it->second->synchronize(data_accessor, tag);
    } else {
      auto synchronizer_it = dummy_synchronizers.find(neighborhood_id);
      if (synchronizer_it == dummy_synchronizers.end())
        continue;

      synchronizer_it->second->synchronizeOnce(data_accessor, tag);
    }
  }
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::averageInternals(const GhostType & ghost_type) {
  /// update the weights of the weight function
  if (ghost_type == _not_ghost)
    this->computeWeights();

  /// loop over all neighborhoods and compute the non-local variables
  for (auto & neighborhood : neighborhoods) {
    /// loop over all the non-local variables of the given neighborhood
    for (auto & non_local_variable : non_local_variables) {
      NonLocalVariable & non_local_var = *non_local_variable.second;
      neighborhood.second->weightedAverageOnNeighbours(
          non_local_var.local, non_local_var.non_local,
          non_local_var.nb_component, ghost_type);
    }
  }
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::computeWeights() {
  AKANTU_DEBUG_IN();

  this->updateWeightFunctionInternals();
  this->volumes.clear();

  for (const auto & global_neighborhood : global_neighborhoods) {
    auto it = neighborhoods.find(global_neighborhood);

    if (it != neighborhoods.end())
      it->second->updateWeights();
    else {
      dummy_synchronizers[global_neighborhood]->synchronize(
          dummy_accessor, SynchronizationTag::_mnl_weight);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::updatePairLists() {
  AKANTU_DEBUG_IN();

  integration_points_positions.initialize(
      this->model.getFEEngine(), _nb_component = spatial_dimension,
      _spatial_dimension = spatial_dimension);

  /// compute the position of the quadrature points
  this->model.getFEEngine().computeIntegrationPointsCoordinates(
      integration_points_positions);

  for (auto & pair : neighborhoods)
    pair.second->updatePairList();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::registerNonLocalVariable(const ID & variable_name,
                                               const ID & nl_variable_name,
                                               UInt nb_component) {

  AKANTU_DEBUG_IN();

  auto non_local_variable_it = non_local_variables.find(variable_name);

  if (non_local_variable_it == non_local_variables.end())
    non_local_variables[nl_variable_name] = std::make_unique<NonLocalVariable>(
        variable_name, nl_variable_name, this->id, nb_component);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
ElementTypeMapReal &
NonLocalManager::registerWeightFunctionInternal(const ID & field_name) {

  AKANTU_DEBUG_IN();

  auto it = this->weight_function_internals.find(field_name);
  if (it == weight_function_internals.end()) {
    weight_function_internals[field_name] =
        std::make_unique<ElementTypeMapReal>(field_name, this->id,
                                             this->memory_id);
  }

  AKANTU_DEBUG_OUT();

  return *(weight_function_internals[field_name]);
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::updateWeightFunctionInternals() {
  for (auto & pair : this->weight_function_internals) {
    auto & internals = *pair.second;
    internals.clear();
    for (auto ghost_type : ghost_types)
      this->callback->updateLocalInternal(internals, ghost_type, _ek_regular);
  }
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::initNonLocalVariables() {
  /// loop over all the non-local variables
  for (auto & pair : non_local_variables) {
    auto & variable = *pair.second;
    variable.non_local.initialize(this->model.getFEEngine(),
                                  _nb_component = variable.nb_component,
                                  _spatial_dimension = spatial_dimension);
  }
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::computeAllNonLocalStresses() {

  /// update the flattened version of the internals
  for (auto & pair : non_local_variables) {
    auto & variable = *pair.second;
    variable.local.clear();
    variable.non_local.clear();
    for (auto ghost_type : ghost_types) {
      this->callback->updateLocalInternal(variable.local, ghost_type,
                                          _ek_regular);
    }
  }

  this->volumes.clear();

  for (auto & pair : neighborhoods) {
    auto & neighborhood = *pair.second;
    neighborhood.asynchronousSynchronize(SynchronizationTag::_mnl_for_average);
  }

  this->averageInternals(_not_ghost);

  AKANTU_DEBUG_INFO("Wait distant non local stresses");

  for (auto & pair : neighborhoods) {
    auto & neighborhood = *pair.second;
    neighborhood.waitEndSynchronize(SynchronizationTag::_mnl_for_average);
  }

  this->averageInternals(_ghost);

  /// copy the results in the materials
  for (auto & pair : non_local_variables) {
    auto & variable = *pair.second;
    for (auto ghost_type : ghost_types) {
      this->callback->updateNonLocalInternal(variable.non_local, ghost_type,
                                             _ek_regular);
    }
  }

  this->callback->computeNonLocalStresses(_not_ghost);

  ++this->compute_stress_calls;
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::cleanupExtraGhostElements() {
  // ElementTypeMap<UInt> & nb_ghost_protected) {

  using ElementSet = std::set<Element>;
  ElementSet relevant_ghost_elements;

  /// loop over all the neighborhoods and get their protected ghosts
  for (auto & pair : neighborhoods) {
    auto & neighborhood = *pair.second;
    ElementSet to_keep_per_neighborhood;

    neighborhood.getRelevantGhostElements(to_keep_per_neighborhood);
    relevant_ghost_elements.insert(to_keep_per_neighborhood.begin(), to_keep_per_neighborhood.end());
  }

    for (auto & pair : neighborhoods) {
    auto & neighborhood = *pair.second;
    neighborhood.cleanupExtraGhostElements(relevant_ghost_elements);
  }

  // /// remove all unneccessary ghosts from the mesh
  // /// Create list of element to remove and new numbering for element to keep
  // Mesh & mesh = this->model.getMesh();
  // ElementSet ghost_to_erase;

  // RemovedElementsEvent remove_elem(mesh);
  // auto & new_numberings = remove_elem.getNewNumbering();
  // Element element;
  // element.ghost_type = _ghost;

  // for (auto & type : mesh.elementTypes(spatial_dimension, _ghost)) {
  //   element.type = type;
  //   UInt nb_ghost_elem = mesh.getNbElement(type, _ghost);
  //   // UInt nb_ghost_elem_protected = 0;
  //   // try {
  //   //   nb_ghost_elem_protected = nb_ghost_protected(type, _ghost);
  //   // } catch (...) {
  //   // }

  //   if (!new_numberings.exists(type, _ghost))
  //     new_numberings.alloc(nb_ghost_elem, 1, type, _ghost);
  //   else
  //     new_numberings(type, _ghost).resize(nb_ghost_elem);

  //   Array<UInt> & new_numbering = new_numberings(type, _ghost);
  //   for (UInt g = 0; g < nb_ghost_elem; ++g) {
  //     element.element = g;
  //     if (element.element >= nb_ghost_elem_protected &&
  //         relevant_ghost_elements.find(element) ==
  //             relevant_ghost_elements.end()) {
  //       remove_elem.getList().push_back(element);
  //       new_numbering(element.element) = UInt(-1);
  //     }
  //   }
  //   /// renumber remaining ghosts
  //   UInt ng = 0;
  //   for (UInt g = 0; g < nb_ghost_elem; ++g) {
  //     if (new_numbering(g) != UInt(-1)) {
  //       new_numbering(g) = ng;
  //       ++ng;
  //     }
  //   }
  // }

  // for (auto & type : mesh.elementTypes(spatial_dimension, _not_ghost)) {
  //   UInt nb_elem = mesh.getNbElement(type, _not_ghost);
  //   if (!new_numberings.exists(type, _not_ghost))
  //     new_numberings.alloc(nb_elem, 1, type, _not_ghost);
  //   Array<UInt> & new_numbering = new_numberings(type, _not_ghost);
  //   for (UInt e = 0; e < nb_elem; ++e) {
  //     new_numbering(e) = e;
  //   }
  // }
  // mesh.sendEvent(remove_elem);
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::onElementsRemoved(
    const Array<Element> & element_list,
    const ElementTypeMapArray<UInt> & new_numbering,
    __attribute__((unused)) const RemovedElementsEvent & event) {

  FEEngine & fee = this->model.getFEEngine();
  this->removeIntegrationPointsFromMap(
      event.getNewNumbering(), spatial_dimension, integration_points_positions,
      fee, _ek_regular);
  this->removeIntegrationPointsFromMap(event.getNewNumbering(), 1, volumes, fee,
                                       _ek_regular);

  /// loop over all the neighborhoods and call onElementsRemoved
  auto global_neighborhood_it = global_neighborhoods.begin();
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
void NonLocalManager::onElementsAdded(const Array<Element> &,
                                      const NewElementsEvent &) {
  this->resizeElementTypeMap(1, volumes, model.getFEEngine());
  this->resizeElementTypeMap(spatial_dimension, integration_points_positions,
                             model.getFEEngine());
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::resizeElementTypeMap(UInt nb_component,
                                           ElementTypeMapReal & element_map,
                                           const FEEngine & fee,
                                           const ElementKind el_kind) {
  Mesh & mesh = this->model.getMesh();

  for (auto gt : ghost_types) {
    for (auto type : mesh.elementTypes(spatial_dimension, gt, el_kind)) {
      UInt nb_element = mesh.getNbElement(type, gt);
      UInt nb_quads = fee.getNbIntegrationPoints(type, gt);
      if (!element_map.exists(type, gt))
        element_map.alloc(nb_element * nb_quads, nb_component, type, gt);
      else
        element_map(type, gt).resize(nb_element * nb_quads);
    }
  }
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::removeIntegrationPointsFromMap(
    const ElementTypeMapArray<UInt> & new_numbering, UInt nb_component,
    ElementTypeMapReal & element_map, const FEEngine & fee,
    const ElementKind el_kind) {

  for (auto gt : ghost_types) {
    for (auto type : new_numbering.elementTypes(_all_dimensions, gt, el_kind)) {
      if (element_map.exists(type, gt)) {
        const Array<UInt> & renumbering = new_numbering(type, gt);

        Array<Real> & vect = element_map(type, gt);
        UInt nb_quad_per_elem = fee.getNbIntegrationPoints(type, gt);
        Array<Real> tmp(renumbering.size() * nb_quad_per_elem, nb_component);

        AKANTU_DEBUG_ASSERT(
            tmp.size() == vect.size(),
            "Something strange append some mater was created or disappeared in "
                << vect.getID() << "(" << vect.size() << "!=" << tmp.size()
                << ") "
                   "!!");

        UInt new_size = 0;
        for (UInt i = 0; i < renumbering.size(); ++i) {
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

/* -------------------------------------------------------------------------- */
UInt NonLocalManager::getNbData(const Array<Element> & elements,
                                const ID & id) const {
  UInt size = 0;
  UInt nb_quadrature_points = this->model.getNbIntegrationPoints(elements);
  auto it = non_local_variables.find(id);

  AKANTU_DEBUG_ASSERT(it != non_local_variables.end(),
                      "The non-local variable " << id << " is not registered");

  size += it->second->nb_component * sizeof(Real) * nb_quadrature_points;
  return size;
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::packData(CommunicationBuffer & buffer,
                               const Array<Element> & elements,
                               const ID & id) const {

  auto it = non_local_variables.find(id);

  AKANTU_DEBUG_ASSERT(it != non_local_variables.end(),
                      "The non-local variable " << id << " is not registered");

  DataAccessor<Element>::packElementalDataHelper<Real>(
      it->second->local, buffer, elements, true, this->model.getFEEngine());
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::unpackData(CommunicationBuffer & buffer,
                                 const Array<Element> & elements,
                                 const ID & id) const {
  auto it = non_local_variables.find(id);

  AKANTU_DEBUG_ASSERT(it != non_local_variables.end(),
                      "The non-local variable " << id << " is not registered");

  DataAccessor<Element>::unpackElementalDataHelper<Real>(
      it->second->local, buffer, elements, true, this->model.getFEEngine());
}

} // namespace akantu
