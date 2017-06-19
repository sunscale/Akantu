/**
 * @file   non_local_neighborhood_base.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Sat Sep 26 2015
 * @date last modification: Wed Nov 25 2015
 *
 * @brief  Implementation of non-local neighborhood base
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
#include "non_local_neighborhood_base.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
NonLocalNeighborhoodBase::NonLocalNeighborhoodBase(const SolidMechanicsModel & model,
						   const ElementTypeMapReal & quad_coordinates,
						   const ID & id,
						   const MemoryID & memory_id)  :
  NeighborhoodBase(model, quad_coordinates, id, memory_id),
  Parsable(_st_non_local, id) {

  AKANTU_DEBUG_IN();

  this->registerParam("radius"       , neighborhood_radius             , 100.,
		      _pat_parsable | _pat_readable  , "Non local radius");

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
NonLocalNeighborhoodBase::~NonLocalNeighborhoodBase() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NonLocalNeighborhoodBase::createGridSynchronizer() {
  this->is_creating_grid = true;
  std::set<SynchronizationTag> tags;
  tags.insert(_gst_mnl_for_average);
  tags.insert(_gst_mnl_weight);

  std::stringstream sstr; sstr << getID() << ":grid_synchronizer";
  this->grid_synchronizer = GridSynchronizer::createGridSynchronizer(this->model.getMesh(),
								     *spatial_grid,
								     sstr.str(),
								     synch_registry,
								     tags, 0, false);
  this->is_creating_grid = false;

}


/* -------------------------------------------------------------------------- */
void NonLocalNeighborhoodBase::cleanupExtraGhostElements(std::set<Element> & relevant_ghost_elements) {

  PairList::const_iterator first_pair = pair_list[_ghost].begin();
  PairList::const_iterator last_pair  = pair_list[_ghost].end();
  for(;first_pair != last_pair; ++first_pair) {
    const IntegrationPoint & q2 = first_pair->second;
    relevant_ghost_elements.insert(q2);
  }

  Array<Element> ghosts_to_erase(0);
  Mesh & mesh = this->model.getMesh();
  Mesh::type_iterator it        = mesh.firstType(spatial_dimension, _ghost);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, _ghost);

  Element element;  
  element.ghost_type = _ghost;

  std::set<Element>::const_iterator end = relevant_ghost_elements.end();
  for(; it != last_type; ++it) {
    element.type = *it;
    UInt nb_ghost_elem = mesh.getNbElement(*it, _ghost);
    for (UInt g = 0; g < nb_ghost_elem; ++g) {
      element.element = g;
      if (relevant_ghost_elements.find(element) == end) {
	ghosts_to_erase.push_back(element);
      }
    }
  }

  /// remove the unneccessary ghosts from the synchronizer
  this->grid_synchronizer->removeElements(ghosts_to_erase);
}

} // akantu
