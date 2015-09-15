/**
 * @file   dof_manager.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Wed Aug 12 09:52:30 2015
 *
 * @brief  Implementation of the common parts of the DOFManagers
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
#include "dof_manager.hh"
#include "mesh.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
DOFManager::DOFManager(const Mesh & mesh, SolverCallback & solver_callback,
                       const ID & id, const MemoryID & memory_id)
  : Memory(id, memory_id), mesh(mesh), solver_callback(solver_callback) {}

/* -------------------------------------------------------------------------- */
void DOFManager::assembleElementalArrayLocalArray(
    const Array<Real> & elementary_vect, Array<Real> & array_assembeled,
    const ElementType & type, const GhostType & ghost_type, Real scale_factor,
    const Array<UInt> & filter_elements) {
  AKANTU_DEBUG_IN();

  UInt nb_element;
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt nb_degree_of_freedom =
      elementary_vect.getNbComponent() / nb_nodes_per_element;

  UInt * filter_it = NULL;
  if (filter_elements != empty_filter) {
    nb_element = filter_elements.getSize();
    filter_it = filter_elements.storage();
  } else {
    nb_element = mesh.getNbElement(type, ghost_type);
  }

  AKANTU_DEBUG_ASSERT(elementary_vect.getSize() == nb_element,
                      "The vector elementary_vect("
                          << elementary_vect.getID()
                          << ") has not the good size.");

  const Array<UInt> connectivity = this->mesh.getConnectivity(type, ghost_type);
  Array<UInt>::const_vector_iterator conn_begin =
      connectivity.begin(nb_nodes_per_element);
  Array<UInt>::const_vector_iterator conn_it = conn_begin;

  UInt size_mat = nb_nodes_per_element * nb_degree_of_freedom;
  Array<Real>::const_vector_iterator elem_it = elementary_vect.begin(size_mat);

  for (UInt el = 0; el < nb_element; ++el, ++elem_it) {
    if (filter_it != NULL)
      conn_it = conn_begin + *filter_it;

    for (UInt n = 0, ld = 0; n < nb_nodes_per_element; ++n) {
      UInt offset_node = (*conn_it)(n) * nb_degree_of_freedom;

      for (UInt d = 0; d < nb_degree_of_freedom; ++d, ++ld) {
        array_assembeled[offset_node + d] += scale_factor * (*elem_it)(ld);
      }
    }

    if (filter_it != NULL)
      ++filter_it;
    else
      ++conn_it;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DOFManager::assembleElementalArrayResidual(
    const ID & dof_id, const Array<Real> & elementary_vect,
    const ElementType & type, const GhostType & ghost_type, Real scale_factor,
    const Array<UInt> & filter_elements) {
  AKANTU_DEBUG_IN();

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt nb_degree_of_freedom =
      elementary_vect.getNbComponent() / nb_nodes_per_element;
  Array<Real> array_localy_assembeled(this->mesh.getNbNodes(),
                                      nb_degree_of_freedom);

  this->assembleElementalArrayLocalArray(
      elementary_vect, array_localy_assembeled, type, ghost_type, scale_factor,
      filter_elements);

  this->assembleToResidual(dof_id, array_localy_assembeled, scale_factor);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DOFManager::registerDOFs(const ID & dof_id, Array<Real> & dofs_array) {
  DOFStorage::iterator it = this->dofs.find(dof_id);

  if (it != this->dofs.end()) {
    AKANTU_EXCEPTION("This dof array has already been registered");
  }

  DOFData * dofs_storage = new DOFData();

  dofs_storage->dof = &dofs_array;
  dofs_storage->blocked_dofs = NULL;

  dofs_storage->local_equation_number.resize(dofs_array.getSize() *
                                             dofs_array.getNbComponent());
  dofs_storage->global_equation_number.resize(dofs_array.getSize() *
                                              dofs_array.getNbComponent());

  this->dofs[dof_id] = dofs_storage;
}

/* -------------------------------------------------------------------------- */
void DOFManager::registerDOFDerivative(const ID & dof_id, UInt order,
                                       Array<Real> & dofs_derivative) {
  DOFStorage::iterator it = this->dofs.find(dof_id);

  if (it == this->dofs.end()) {
    AKANTU_EXCEPTION("The dof array corresponding to this derivatives has not "
                     << "been registered yet");
  }

  DOFData & dof = *it->second;
  std::vector<Array<Real> *> & derivatives = dof.dof_derivatives;

  if (derivatives.size() < order) {
    derivatives.resize(order, NULL);
  } else {
    if (derivatives[order - 1] != NULL) {
      AKANTU_EXCEPTION("The dof derivatives of order "
                       << order << " already been registered for this dof ("
                       << dof_id << ")");
    }
  }

  derivatives[order - 1] = &dofs_derivative;
}

/* -------------------------------------------------------------------------- */
void DOFManager::registerBlockedDOFs(const ID & dof_id,
                                     Array<Real> & blocked_dofs) {
  DOFStorage::iterator it = this->dofs.find(dof_id);

  if (it == this->dofs.end()) {
    AKANTU_EXCEPTION("The dof array corresponding to this derivatives has not "
                     << "been registered yet");
  }

  DOFData & dof = *it->second;

  if (dof.blocked_dofs != NULL) {
    AKANTU_EXCEPTION("The blocked dofs array for "
                     << dof_id << " has already been registered");
  }

  dof.blocked_dofs = &blocked_dofs;
}

/* -------------------------------------------------------------------------- */

__END_AKANTU__
