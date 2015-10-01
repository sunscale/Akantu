/**
 * @file   non_local_manager.cc
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Mon Sep 21 15:32:10 2015
 *
 * @brief  Implementation of non-local manager
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
#include "non_local_manager.hh"
#include "non_local_neighborhood.hh"
/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
NonLocalManager::NonLocalManager(SolidMechanicsModel & model, 
				const ID & id,
				const MemoryID & memory_id) : 
  Memory(id, memory_id),
  model(model),
  volumes("volumes", id){ 
  UInt spatial_dimension = this->model.getSpatialDimension();
  Mesh & mesh = this->model.getMesh();

  mesh.initElementTypeMapArray(volumes, 1, spatial_dimension, false, _ek_regular, true);

  for(UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType) g;
    Mesh::type_iterator it  = mesh.firstType(spatial_dimension, gt, _ek_regular);
    Mesh::type_iterator end = mesh.lastType(spatial_dimension, gt, _ek_regular);
    for(; it != end; ++it) {
      UInt nb_element = mesh.getNbElement(*it, gt);
      UInt nb_quads = this->model.getFEEngine().getNbQuadraturePoints(*it, gt);
      volumes.alloc(nb_element *nb_quads, 1, *it, gt);
    }
  }
#ifdef AKANTU_IGFEM
  mesh.initElementTypeMapArray(volumes, 1, spatial_dimension, false, _ek_igfem, true);
  for(UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType) g;
    Mesh::type_iterator it  = mesh.firstType(spatial_dimension, gt, _ek_igfem);
    Mesh::type_iterator end = mesh.lastType(spatial_dimension, gt, _ek_igfem);
    for(; it != end; ++it) {
      UInt nb_element = mesh.getNbElement(*it, gt);
      UInt nb_quads = this->model.getFEEngine("IGFEMFEEngine").getNbQuadraturePoints(*it, gt);
      volumes.alloc(nb_element *nb_quads, 1, *it, gt);
    }
  }

#endif

}

/* -------------------------------------------------------------------------- */
NonLocalManager::~NonLocalManager() {
  NeighborhoodMap::iterator it;
  for (it = neighborhoods.begin(); it != neighborhoods.end(); ++it) {
    if(it->second) delete it->second;
  }

}

/* -------------------------------------------------------------------------- */
void NonLocalManager::setJacobians(const FEEngine & fe_engine, const ElementKind & kind) {
  Mesh & mesh = this->model.getMesh();
  UInt spatial_dimension = this->model.getSpatialDimension();
  for(UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType) g;
    Mesh::type_iterator it = mesh.firstType(spatial_dimension, gt, kind);
    Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, gt, kind);
    for(; it != last_type; ++it) {
      jacobians(*it, gt) = &fe_engine.getIntegratorInterface().getJacobians(*it, gt);
    }
  }
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::createNeighborhood(const ID & type, Real radius, const ID & name) {

  AKANTU_DEBUG_IN();

  /// check if neighborhood already exists
  NeighborhoodMap::const_iterator it = neighborhoods.find(name);

  if (it == neighborhoods.end() ) {
    /// create new neighborhood for given ID
    std::stringstream sstr; sstr << "neighborhood:" << name;
    neighborhoods[name] = new NonLocalNeighborhood(*this, radius, type, sstr.str());
  }

  AKANTU_DEBUG_OUT();

}

/* -------------------------------------------------------------------------- */
void NonLocalManager::createNeighborhoodSynchronizers() {

  NeighborhoodMap::iterator it;
  for (it = neighborhoods.begin(); it != neighborhoods.end(); ++it) {
    it->second->createGridSynchronizer();
  }

}


__END_AKANTU__
