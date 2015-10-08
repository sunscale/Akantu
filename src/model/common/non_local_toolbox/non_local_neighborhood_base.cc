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
								     tags);
  this->is_creating_grid = false;

}


__END_AKANTU__
