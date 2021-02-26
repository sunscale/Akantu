/**
 * @file   non_local_manager_igfem.cc
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Mon Sep 21 15:32:10 2015
 *
 * @brief  Implementation of non-local manager igfem
 *
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
#ifdef AKANTU_DAMAGE_NON_LOCAL
#include "non_local_manager_igfem.hh"
#include "material_non_local.hh"
/* -------------------------------------------------------------------------- */
namespace akantu {

/* -------------------------------------------------------------------------- */
NonLocalManagerIGFEM::NonLocalManagerIGFEM(SolidMechanicsModelIGFEM & model,
                                           const ID & id,
                                           const MemoryID & memory_id)
    : NonLocalManager(model, id, memory_id) {

  Mesh & mesh = this->model.getMesh();

  /// initialize the element type map array
  /// it will be resized to nb_quad * nb_element during the computation of
  /// coords
  mesh.initElementTypeMapArray(quad_positions, spatial_dimension,
                               spatial_dimension, false, _ek_igfem, true);
}

/* -------------------------------------------------------------------------- */
NonLocalManagerIGFEM::~NonLocalManagerIGFEM() {}

/* -------------------------------------------------------------------------- */
void NonLocalManagerIGFEM::init() {

  /// store the number of current ghost elements for each type in the mesh
  ElementTypeMap<UInt> nb_ghost_protected;
  Mesh & mesh = this->model.getMesh();
  for (UInt k = _ek_regular; k <= _ek_igfem; ++k) {
    ElementKind el_kind = (ElementKind)k;
    Mesh::type_iterator it = mesh.firstType(spatial_dimension, _ghost, el_kind);
    Mesh::type_iterator last_type =
        mesh.lastType(spatial_dimension, _ghost, el_kind);
    for (; it != last_type; ++it)
      nb_ghost_protected(mesh.getNbElement(*it, _ghost), *it, _ghost);
  }

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

  FEEngine & fee_regular = this->model.getFEEngine();
  FEEngine & fee_igfem = this->model.getFEEngine("IGFEMFEEngine");

  this->updatePairLists();
  /// cleanup the unneccessary ghost elements
  this->cleanupExtraGhostElements(nb_ghost_protected);
  this->initElementTypeMap(1, volumes, fee_regular, _ek_regular);
  this->initElementTypeMap(1, volumes, fee_igfem, _ek_igfem);
  this->setJacobians(fee_regular, _ek_regular);
  this->setJacobians(fee_igfem, _ek_igfem);
  this->initNonLocalVariables();
  this->computeWeights();
}

/* -------------------------------------------------------------------------- */
void NonLocalManagerIGFEM::computeAllNonLocalStresses() {

  /// update the flattened version of the internals
  std::map<ID, NonLocalVariable *>::iterator non_local_variable_it =
      non_local_variables.begin();
  std::map<ID, NonLocalVariable *>::iterator non_local_variable_end =
      non_local_variables.end();

  for (; non_local_variable_it != non_local_variable_end;
       ++non_local_variable_it) {
    non_local_variable_it->second->local.zero();
    non_local_variable_it->second->non_local.zero();
    for (UInt gt = _not_ghost; gt <= _ghost; ++gt) {
      GhostType ghost_type = (GhostType)gt;
      this->flattenInternal(non_local_variable_it->second->local, ghost_type,
                            _ek_regular);
      this->flattenInternal(non_local_variable_it->second->local, ghost_type,
                            _ek_igfem);
    }
  }

  this->volumes.zero();
  /// loop over all the neighborhoods and compute intiate the
  /// exchange of the non-local_variables
  std::set<ID>::const_iterator global_neighborhood_it =
      global_neighborhoods.begin();
  NeighborhoodMap::iterator it;
  for (; global_neighborhood_it != global_neighborhoods.end();
       ++global_neighborhood_it) {
    it = neighborhoods.find(*global_neighborhood_it);
    if (it != neighborhoods.end())
      it->second->getSynchronizerRegistry().asynchronousSynchronize(
          SynchronizationTag::_mnl_for_average);
    else
      dummy_synchronizers[*global_neighborhood_it]->asynchronousSynchronize(
          dummy_accessor, SynchronizationTag::_mnl_for_average);
  }

  this->averageInternals(_not_ghost);

  AKANTU_DEBUG_INFO("Wait distant non local stresses");

  /// loop over all the neighborhoods and block until all non-local
  /// variables have been exchanged
  global_neighborhood_it = global_neighborhoods.begin();
  it = neighborhoods.begin();
  for (; global_neighborhood_it != global_neighborhoods.end();
       ++global_neighborhood_it) {
    it = neighborhoods.find(*global_neighborhood_it);
    if (it != neighborhoods.end())
      it->second->getSynchronizerRegistry().waitEndSynchronize(
          SynchronizationTag::_mnl_for_average);
    else
      dummy_synchronizers[*global_neighborhood_it]->waitEndSynchronize(
          dummy_accessor, SynchronizationTag::_mnl_for_average);
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
void NonLocalManagerIGFEM::cleanupExtraGhostElements(
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
    to_keep_per_neighborhood.zero();
  }

  /// remove all unneccessary ghosts from the mesh
  /// Create list of element to remove and new numbering for element to keep
  Mesh & mesh = this->model.getMesh();
  ElementSet ghost_to_erase;
  RemovedElementsEvent remove_elem(mesh);
  Element element;

  for (UInt k = _ek_regular; k < _ek_igfem; ++k) {
    ElementKind el_kind = (ElementKind)k;
    element.kind = _ek_igfem;
    Mesh::type_iterator it = mesh.firstType(spatial_dimension, _ghost, el_kind);
    Mesh::type_iterator last_type =
        mesh.lastType(spatial_dimension, _ghost, el_kind);

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
  }

  for (UInt k = _ek_regular; k < _ek_igfem; ++k) {
    ElementKind el_kind = (ElementKind)k;
    Mesh::type_iterator it = mesh.firstType(spatial_dimension, _ghost, el_kind);
    Mesh::type_iterator last_type =
        mesh.lastType(spatial_dimension, _ghost, el_kind);
    for (; it != last_type; ++it) {
      UInt nb_elem = mesh.getNbElement(*it, _not_ghost);
      if (!remove_elem.getNewNumbering().exists(*it, _not_ghost))
        remove_elem.getNewNumbering().alloc(nb_elem, 1, *it, _not_ghost);
      Array<UInt> & new_numbering =
          remove_elem.getNewNumbering(*it, _not_ghost);
      for (UInt e = 0; e < nb_elem; ++e) {
        new_numbering(e) = e;
      }
    }
  }

  mesh.sendEvent(remove_elem);
}

/* -------------------------------------------------------------------------- */
void NonLocalManagerIGFEM::onElementsAdded(__attribute__((unused))
                                           const Array<Element> & element_list,
                                           __attribute__((unused))
                                           const NewElementsEvent & event) {

  FEEngine & fee = this->model.getFEEngine("IGFEMFEEngine");
  this->resizeElementTypeMap(1, volumes, fee, _ek_igfem);
  this->resizeElementTypeMap(spatial_dimension, quad_positions, fee, _ek_igfem);

  NonLocalManager::onElementsAdded(element_list, event);
}

/* -------------------------------------------------------------------------- */
void NonLocalManagerIGFEM::onElementsRemoved(
    const Array<Element> & element_list,
    const ElementTypeMapArray<UInt> & new_numbering,
    __attribute__((unused)) const RemovedElementsEvent & event) {

  FEEngine & fee = this->model.getFEEngine("IGFEMFEEngine");

  this->removeIntegrationPointsFromMap(event.getNewNumbering(),
                                       spatial_dimension, quad_positions, fee,
                                       _ek_igfem);
  this->removeIntegrationPointsFromMap(event.getNewNumbering(), 1, volumes, fee,
                                       _ek_igfem);

  NonLocalManager::onElementsRemoved(element_list, new_numbering, event);
}

} // namespace akantu

#endif /* AKANTU_DAMAGE_NON_LOCAL */
