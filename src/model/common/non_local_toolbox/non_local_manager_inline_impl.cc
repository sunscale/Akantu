/**
 * @file   non_local_manager_inline_impl.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Wed Nov 25 2015
 *
 * @brief  inline implementation of non-local manager functions
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
#include "grid_synchronizer.hh"
#include "synchronizer_registry.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_NON_LOCAL_MANAGER_INLINE_IMPL_CC__
#define __AKANTU_NON_LOCAL_MANAGER_INLINE_IMPL_CC__

namespace akantu {

/* -------------------------------------------------------------------------- */
inline void NonLocalManager::insertQuad(const IntegrationPoint & quad,
                                        const Vector<Real> & coords,
                                        const ID & neighborhood) {

  AKANTU_DEBUG_ASSERT(neighborhoods[neighborhood] != NULL,
                      "The neighborhood " << neighborhood
                                          << "has not been created");

  neighborhoods[neighborhood]->insertQuad(quad, coords);
}

/* -------------------------------------------------------------------------- */
inline void NonLocalManager::registerNeighborhood(const ID & neighborhood,
                                                  const ID & weight_func_id) {

  /// check if neighborhood has already been created
  NeighborhoodMap::const_iterator it = neighborhoods.find(neighborhood);
  if (it == neighborhoods.end()) {
    this->createNeighborhood(weight_func_id, neighborhood);
  }
}

/* -------------------------------------------------------------------------- */
inline NonLocalNeighborhoodBase &
NonLocalManager::getNeighborhood(const ID & name) const {
  AKANTU_DEBUG_IN();

  NeighborhoodMap::const_iterator it = neighborhoods.find(name);

  AKANTU_DEBUG_ASSERT(it != neighborhoods.end(),
                      "The neighborhood " << name << " is not registered");

  AKANTU_DEBUG_OUT();
  return *(it->second);
}

/* -------------------------------------------------------------------------- */
inline void NonLocalManager::computeWeights() {
  AKANTU_DEBUG_IN();

  this->updateWeightFunctionInternals();
  this->volumes.clear();

  // NeighborhoodMap::iterator it = neighborhoods.begin();
  // NeighborhoodMap::iterator end = neighborhoods.end();

  // for (; it != end; ++it)
  //   it->second->updateWeights();

  /// loop over all the neighborhoods and call onElementsRemoved
  std::set<ID>::const_iterator global_neighborhood_it =
      global_neighborhoods.begin();
  NeighborhoodMap::iterator it;
  for (; global_neighborhood_it != global_neighborhoods.end();
       ++global_neighborhood_it) {
    it = neighborhoods.find(*global_neighborhood_it);
    if (it != neighborhoods.end())
      it->second->updateWeights();
    else {
      dummy_synchronizers[*global_neighborhood_it]->asynchronousSynchronize(
          dummy_accessor, _gst_mnl_weight);
      dummy_synchronizers[*global_neighborhood_it]->waitEndSynchronize(
          dummy_accessor, _gst_mnl_weight);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline void NonLocalManager::updatePairLists() {
  AKANTU_DEBUG_IN();

  /// compute the position of the quadrature points
  this->model.getFEEngine().computeIntegrationPointsCoordinates(quad_positions);

  NeighborhoodMap::iterator it = neighborhoods.begin();
  NeighborhoodMap::iterator end = neighborhoods.end();

  for (; it != end; ++it)
    it->second->updatePairList();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline void NonLocalManager::registerNonLocalVariable(
    const ID & variable_name, const ID & nl_variable_name, UInt nb_component) {

  AKANTU_DEBUG_IN();

  std::map<ID, NonLocalVariable *>::iterator non_local_variable_it =
      non_local_variables.find(variable_name);

  if (non_local_variable_it == non_local_variables.end())
    non_local_variables[nl_variable_name] = new NonLocalVariable(
        variable_name, nl_variable_name, this->id, nb_component);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline void NonLocalManager::registerNonLocalMaterial(Material & new_mat) {
  non_local_materials.push_back(&new_mat);
}

/* -------------------------------------------------------------------------- */
inline ElementTypeMapReal &
NonLocalManager::registerWeightFunctionInternal(const ID & field_name) {

  AKANTU_DEBUG_IN();

  std::map<ID, ElementTypeMapReal *>::const_iterator it =
      this->weight_function_internals.find(field_name);
  if (it == weight_function_internals.end()) {
    weight_function_internals[field_name] = new ElementTypeMapReal(field_name, id, memory_id);
  }

  return *(weight_function_internals[field_name]);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline void NonLocalManager::updateWeightFunctionInternals() {

  std::map<ID, ElementTypeMapReal *>::const_iterator it =
      this->weight_function_internals.begin();
  std::map<ID, ElementTypeMapReal *>::const_iterator end =
      this->weight_function_internals.end();
  for (; it != end; ++it) {
    it->second->clear();
    for (UInt g = _not_ghost; g <= _ghost; ++g) {
      GhostType ghost_type = (GhostType)g;
      this->flattenInternal(*(it->second), ghost_type, _ek_regular);
    }
  }
}

/* -------------------------------------------------------------------------- */
inline void
NonLocalManager::nonLocalVariableToNeighborhood(const ID & variable_name,
                                                const ID & neighborhood) {

  NeighborhoodMap::const_iterator it = neighborhoods.find(neighborhood);

  AKANTU_DEBUG_ASSERT(it != neighborhoods.end(), "The neighborhood "
                                                     << neighborhood
                                                     << " is not registered");
  it->second->registerNonLocalVariable(variable_name);
}

/* -------------------------------------------------------------------------- */
inline UInt
NonLocalManager::getNbDataForElements(const Array<Element> & elements,
                                      const ID & id) const {
  UInt size = 0;
  UInt nb_quadrature_points = this->getModel().getNbIntegrationPoints(elements);
  std::map<ID, NonLocalVariable *>::const_iterator it =
      non_local_variables.find(id);

  AKANTU_DEBUG_ASSERT(it != non_local_variables.end(),
                      "The non-local variable " << id << " is not registered");

  size += it->second->nb_component * sizeof(Real) * nb_quadrature_points;
  return size;
}

/* -------------------------------------------------------------------------- */
inline void NonLocalManager::packElementData(CommunicationBuffer & buffer,
                                             const Array<Element> & elements,
                                             __attribute__((unused))
                                             SynchronizationTag tag,
                                             const ID & id) const {

  std::map<ID, NonLocalVariable *>::const_iterator it =
      non_local_variables.find(id);

  AKANTU_DEBUG_ASSERT(it != non_local_variables.end(),
                      "The non-local variable " << id << " is not registered");

  DataAccessor::packElementalDataHelper<Real>(
      it->second->local, buffer, elements, true, this->model.getFEEngine());
}

/* -------------------------------------------------------------------------- */
inline void NonLocalManager::unpackElementData(CommunicationBuffer & buffer,
                                               const Array<Element> & elements,
                                               __attribute__((unused))
                                               SynchronizationTag tag,
                                               const ID & id) const {

  std::map<ID, NonLocalVariable *>::const_iterator it =
      non_local_variables.find(id);

  AKANTU_DEBUG_ASSERT(it != non_local_variables.end(),
                      "The non-local variable " << id << " is not registered");

  DataAccessor::unpackElementalDataHelper<Real>(
      it->second->local, buffer, elements, true, this->model.getFEEngine());
}

} // akantu

#endif /* __AKANTU_NON_LOCAL_MANAGER_INLINE_IMPL_CC__ */
