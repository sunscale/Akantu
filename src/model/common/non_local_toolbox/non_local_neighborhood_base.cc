/**
 * @file   non_local_neighborhood_base.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Sat Sep 26 2015
 * @date last modification: Fri Dec 08 2017
 *
 * @brief  Implementation of non-local neighborhood base
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
#include "non_local_neighborhood_base.hh"
#include "grid_synchronizer.hh"
#include "model.hh"
/* -------------------------------------------------------------------------- */
#include <memory>

namespace akantu {

/* -------------------------------------------------------------------------- */
NonLocalNeighborhoodBase::NonLocalNeighborhoodBase(
    Model & model, const ElementTypeMapReal & quad_coordinates, const ID & id,
    const MemoryID & memory_id)
    : NeighborhoodBase(model, quad_coordinates, id, memory_id),
      Parsable(ParserType::_non_local, id) {

  AKANTU_DEBUG_IN();

  this->registerParam("radius", neighborhood_radius, 100.,
                      _pat_parsable | _pat_readable, "Non local radius");

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
NonLocalNeighborhoodBase::~NonLocalNeighborhoodBase() = default;

/* -------------------------------------------------------------------------- */
void NonLocalNeighborhoodBase::createGridSynchronizer() {
  this->is_creating_grid = true;

  this->grid_synchronizer = std::make_unique<GridSynchronizer>(
      this->model.getMesh(), *spatial_grid, *this,
      std::set<SynchronizationTag>{SynchronizationTag::_mnl_weight,
                                   SynchronizationTag::_mnl_for_average},
      std::string(getID() + ":grid_synchronizer"), this->memory_id, false);

  this->is_creating_grid = false;
}

/* -------------------------------------------------------------------------- */
void NonLocalNeighborhoodBase::synchronize(
    DataAccessor<Element> & data_accessor, const SynchronizationTag & tag) {
  if (not grid_synchronizer) {
    return;
  }

  grid_synchronizer->synchronizeOnce(data_accessor, tag);
}

/* -------------------------------------------------------------------------- */
void NonLocalNeighborhoodBase::getRelevantGhostElements(
    std::set<Element> & relevant_ghost_elements) {

  for (auto && ghost_type : ghost_type_t{}) {
    auto & pair_list = this->pair_list.at(ghost_type);
    for (auto && pair : pair_list) {
      if (pair.first.ghost_type == _ghost) {
        relevant_ghost_elements.insert(pair.first);
      }
      if (pair.second.ghost_type == _ghost) {
        relevant_ghost_elements.insert(pair.second);
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
void NonLocalNeighborhoodBase::cleanupExtraGhostElements(
    std::set<Element> & relevant_ghost_elements) {

  Array<Element> ghosts_to_erase;
  auto & mesh = this->model.getMesh();

  auto end = relevant_ghost_elements.end();
  for (const auto & type : mesh.elementTypes(
           _spatial_dimension = spatial_dimension, _ghost_type = _ghost)) {
    auto nb_ghost_elem = mesh.getNbElement(type, _ghost);
    for (UInt g = 0; g < nb_ghost_elem; ++g) {
      Element element{type, g, _ghost};
      if (relevant_ghost_elements.find(element) == end) {
        ghosts_to_erase.push_back(element);
      }
    }
  }

  /// remove the unneccessary ghosts from the synchronizer
  mesh.eraseElements(ghosts_to_erase);
}

/* -------------------------------------------------------------------------- */
void NonLocalNeighborhoodBase::registerNonLocalVariable(const ID & id) {
  this->non_local_variables.insert(id);
}

} // namespace akantu
