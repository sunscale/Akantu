/**
 * @file   non_local_manager_inline_impl.cc
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Mon Sep 21 17:30:03 2015
 *
 * @brief  inline implementation of non-local manager functions
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
#include "grid_synchronizer.hh"
#include "synchronizer_registry.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_NON_LOCAL_MANAGER_INLINE_IMPL_CC__
#define __AKANTU_NON_LOCAL_MANAGER_INLINE_IMPL_CC__

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
inline void NonLocalManager::insertQuad(const QuadraturePoint & quad, const Vector<Real> & coords, 
					Real radius, const ID & type, ID name) {

  if (name == "") name = default_neighborhood;

  NeighborhoodMap::const_iterator it = neighborhoods.find(name);
  if (it == neighborhoods.end()) {
    this->createNeighborhood(type, radius, name);
  }

  neighborhoods[name]->insertQuad(quad, coords);

}

/* -------------------------------------------------------------------------- */
inline NonLocalNeighborhoodBase & NonLocalManager::getNeighborhood(const ID & name) const{
  AKANTU_DEBUG_IN();
  ID tmp_name = name;
  if (name == "") tmp_name = default_neighborhood;

  NeighborhoodMap::const_iterator it = neighborhoods.find(tmp_name);

  AKANTU_DEBUG_ASSERT(it != neighborhoods.end(),
		      "The neighborhood " << tmp_name << " is not registered");

  AKANTU_DEBUG_OUT();
  return *(it->second);
}

/* -------------------------------------------------------------------------- */
inline void NonLocalManager::computeWeights() {
  AKANTU_DEBUG_IN();

  Mesh & mesh = this->model.getMesh();
  UInt spatial_dimension = this->model.getSpatialDimension();
  /// need to resize the arrays
  for(UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType) g;
    Mesh::type_iterator it  = mesh.firstType(spatial_dimension, gt, _ek_regular);
    Mesh::type_iterator end = mesh.lastType(spatial_dimension, gt, _ek_regular);
    for(; it != end; ++it) {
      UInt nb_element = mesh.getNbElement(*it, gt);
      UInt nb_quads = this->model.getFEEngine().getNbQuadraturePoints(*it, gt);
      volumes(*it, gt).resize(nb_element *nb_quads);
    }
  }

  this->volumes.clear();

  NeighborhoodMap::iterator it = neighborhoods.begin();
  NeighborhoodMap::iterator end = neighborhoods.end();

  for (; it != end; ++it)
    it->second->computeWeights();

  AKANTU_DEBUG_OUT();

}

/* -------------------------------------------------------------------------- */
inline void NonLocalManager::updatePairLists() {
  AKANTU_DEBUG_IN();

  /// compute the position of the quadrature points
  this->model.getFEEngine().computeQuadraturePointsCoordinates(quad_positions);

  NeighborhoodMap::iterator it = neighborhoods.begin();
  NeighborhoodMap::iterator end = neighborhoods.end();

  for (; it != end; ++it)
    it->second->updatePairList();

  AKANTU_DEBUG_OUT();

}

/* -------------------------------------------------------------------------- */
inline void NonLocalManager::registerNonLocalVariable(const ID & variable_name, const ID & nl_variable_name, UInt nb_component) {

  AKANTU_DEBUG_IN();

  std::map<ID, NonLocalVariable *>::iterator non_local_variable_it = non_local_variables.find(variable_name);
  
  if (non_local_variable_it == non_local_variables.end()) 
    non_local_variables[variable_name] = new NonLocalVariable(variable_name, nl_variable_name, this->id, nb_component);
  

  AKANTU_DEBUG_OUT();

}

/* -------------------------------------------------------------------------- */
inline void NonLocalManager::registerNonLocalMaterial(Material & new_mat) {
  non_local_materials.push_back(&new_mat);
}


__END_AKANTU__

#endif /* __AKANTU_NON_LOCAL_MANAGER_INLINE_IMPL_CC__ */







