/**
 * @file   dof_manager_inline_impl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Feb 21 2013
 * @date last modification: Wed Jan 31 2018
 *
 * @brief  inline functions of the dof manager
 *
 *
 * Copyright (©) 2014-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "dof_manager.hh"
#include "element_group.hh"
#include "solver_vector.hh"
#include "sparse_matrix.hh"
#include "terms_to_assemble.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_DOF_MANAGER_INLINE_IMPL_HH_
#define AKANTU_DOF_MANAGER_INLINE_IMPL_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
inline bool DOFManager::hasDOFs(const ID & dof_id) const {
  auto it = this->dofs.find(dof_id);
  return it != this->dofs.end();
}

/* -------------------------------------------------------------------------- */
inline DOFManager::DOFData & DOFManager::getDOFData(const ID & dof_id) {
  auto it = this->dofs.find(dof_id);
  if (it == this->dofs.end()) {
    AKANTU_EXCEPTION("The dof " << dof_id << " does not exists in "
                                << this->id);
  }
  return *it->second;
}

/* -------------------------------------------------------------------------- */
const DOFManager::DOFData & DOFManager::getDOFData(const ID & dof_id) const {
  auto it = this->dofs.find(dof_id);
  if (it == this->dofs.end()) {
    AKANTU_EXCEPTION("The dof " << dof_id << " does not exists in "
                                << this->id);
  }
  return *it->second;
}

/* -------------------------------------------------------------------------- */
inline void DOFManager::extractElementEquationNumber(
    const Array<Int> & equation_numbers, const Vector<UInt> & connectivity,
    UInt nb_degree_of_freedom, Vector<Int> & element_equation_number) {
  for (UInt i = 0, ld = 0; i < connectivity.size(); ++i) {
    UInt n = connectivity(i);
    for (UInt d = 0; d < nb_degree_of_freedom; ++d, ++ld) {
      element_equation_number(ld) =
          equation_numbers(n * nb_degree_of_freedom + d);
    }
  }
}

/* -------------------------------------------------------------------------- */
template <class DOFData_>
inline DOFData_ & DOFManager::getDOFDataTyped(const ID & dof_id) {
  return aka::as_type<DOFData_>(this->getDOFData(dof_id));
}

/* -------------------------------------------------------------------------- */
template <class DOFData_>
inline const DOFData_ & DOFManager::getDOFDataTyped(const ID & dof_id) const {
  return aka::as_type<DOFData_>(this->getDOFData(dof_id));
}

/* -------------------------------------------------------------------------- */
inline Array<Real> & DOFManager::getDOFs(const ID & dofs_id) {
  return *(this->getDOFData(dofs_id).dof);
}

/* -------------------------------------------------------------------------- */
inline DOFSupportType DOFManager::getSupportType(const ID & dofs_id) const {
  return this->getDOFData(dofs_id).support_type;
}

/* -------------------------------------------------------------------------- */
inline Array<Real> & DOFManager::getPreviousDOFs(const ID & dofs_id) {
  return *(this->getDOFData(dofs_id).previous);
}

/* -------------------------------------------------------------------------- */
inline bool DOFManager::hasPreviousDOFs(const ID & dofs_id) const {
  return (this->getDOFData(dofs_id).previous != nullptr);
}

/* -------------------------------------------------------------------------- */
inline Array<Real> & DOFManager::getDOFsIncrement(const ID & dofs_id) {
  return *(this->getDOFData(dofs_id).increment);
}

/* -------------------------------------------------------------------------- */
inline bool DOFManager::hasDOFsIncrement(const ID & dofs_id) const {
  return (this->getDOFData(dofs_id).increment != nullptr);
}

/* -------------------------------------------------------------------------- */
inline Array<Real> & DOFManager::getDOFsDerivatives(const ID & dofs_id,
                                                    UInt order) {

  if (order == 0) {
    return getDOFs(dofs_id);
  }

  std::vector<Array<Real> *> & derivatives =
      this->getDOFData(dofs_id).dof_derivatives;
  if ((order > derivatives.size()) || (derivatives[order - 1] == nullptr)) {
    AKANTU_EXCEPTION("No derivatives of order " << order << " present in "
                                                << this->id << " for dof "
                                                << dofs_id);
  }

  return *derivatives[order - 1];
}

/* -------------------------------------------------------------------------- */
inline bool DOFManager::hasDOFsDerivatives(const ID & dofs_id,
                                           UInt order) const {
  const std::vector<Array<Real> *> & derivatives =
      this->getDOFData(dofs_id).dof_derivatives;
  return ((order < derivatives.size()) && (derivatives[order - 1] != nullptr));
}

/* -------------------------------------------------------------------------- */
inline const Array<Real> & DOFManager::getSolution(const ID & dofs_id) const {
  return this->getDOFData(dofs_id).solution;
}

/* -------------------------------------------------------------------------- */
inline Array<Real> & DOFManager::getSolution(const ID & dofs_id) {
  return this->getDOFData(dofs_id).solution;
}

/* -------------------------------------------------------------------------- */
inline const Array<bool> &
DOFManager::getBlockedDOFs(const ID & dofs_id) const {
  return *(this->getDOFData(dofs_id).blocked_dofs);
}

/* -------------------------------------------------------------------------- */
inline bool DOFManager::hasBlockedDOFs(const ID & dofs_id) const {
  return (this->getDOFData(dofs_id).blocked_dofs != nullptr);
}

/* -------------------------------------------------------------------------- */
inline bool DOFManager::isLocalOrMasterDOF(UInt dof_num) {
  auto dof_flag = this->dofs_flag(dof_num);
  return (dof_flag & NodeFlag::_local_master_mask) == NodeFlag::_normal;
}

/* -------------------------------------------------------------------------- */
inline bool DOFManager::isSlaveDOF(UInt dof_num) {
  auto dof_flag = this->dofs_flag(dof_num);
  return (dof_flag & NodeFlag::_shared_mask) == NodeFlag::_slave;
}

/* -------------------------------------------------------------------------- */
inline bool DOFManager::isPureGhostDOF(UInt dof_num) {
  auto dof_flag = this->dofs_flag(dof_num);
  return (dof_flag & NodeFlag::_shared_mask) == NodeFlag::_pure_ghost;
}

/* -------------------------------------------------------------------------- */
inline Int DOFManager::localToGlobalEquationNumber(Int local) const {
  return this->global_equation_number(local);
}

/* -------------------------------------------------------------------------- */
inline bool DOFManager::hasGlobalEquationNumber(Int global) const {
  auto it = this->global_to_local_mapping.find(global);
  return (it != this->global_to_local_mapping.end());
}

/* -------------------------------------------------------------------------- */
inline Int DOFManager::globalToLocalEquationNumber(Int global) const {
  auto it = this->global_to_local_mapping.find(global);
  AKANTU_DEBUG_ASSERT(it != this->global_to_local_mapping.end(),
                      "This global equation number "
                          << global << " does not exists in " << this->id);

  return it->second;
}

/* -------------------------------------------------------------------------- */
inline NodeFlag DOFManager::getDOFFlag(Int local_id) const {
  return this->dofs_flag(local_id);
}

/* -------------------------------------------------------------------------- */
inline const Array<UInt> &
DOFManager::getDOFsAssociatedNodes(const ID & dof_id) const {
  const auto & dof_data = this->getDOFData(dof_id);
  return dof_data.associated_nodes;
}

/* -------------------------------------------------------------------------- */
const Array<Int> &
DOFManager::getLocalEquationsNumbers(const ID & dof_id) const {
  return getDOFData(dof_id).local_equation_number;
}

/* -------------------------------------------------------------------------- */
template <typename Vec>
void DOFManager::assembleMatMulVectToArray_(const ID & dof_id, const ID & A_id,
                                            const Array<Real> & x,
                                            Array<Real> & array,
                                            Real scale_factor) {
  Vec tmp_array(aka::as_type<Vec>(*data_cache), this->id + ":tmp_array");
  tmp_array.zero();
  assembleMatMulVectToGlobalArray(dof_id, A_id, x, tmp_array, scale_factor);
  getArrayPerDOFs(dof_id, tmp_array, array);
}

/* -------------------------------------------------------------------------- */
template <typename Mat>
void DOFManager::assembleElementalMatricesToMatrix_(
    Mat & A, const ID & dof_id, const Array<Real> & elementary_mat,
    ElementType type, GhostType ghost_type,
    const MatrixType & elemental_matrix_type,
    const Array<UInt> & filter_elements) {
  AKANTU_DEBUG_IN();

  auto & dof_data = this->getDOFData(dof_id);

  AKANTU_DEBUG_ASSERT(dof_data.support_type == _dst_nodal,
                      "This function applies only on Nodal dofs");

  const auto & equation_number = this->getLocalEquationsNumbers(dof_id);

  UInt nb_element;
  UInt * filter_it = nullptr;
  if (filter_elements != empty_filter) {
    nb_element = filter_elements.size();
    filter_it = filter_elements.storage();
  } else {
    if (dof_data.group_support != "__mesh__") {
      const auto & group_elements =
          this->mesh->getElementGroup(dof_data.group_support)
              .getElements(type, ghost_type);
      nb_element = group_elements.size();
      filter_it = group_elements.storage();
    } else {
      nb_element = this->mesh->getNbElement(type, ghost_type);
    }
  }

  AKANTU_DEBUG_ASSERT(elementary_mat.size() == nb_element,
                      "The vector elementary_mat("
                          << elementary_mat.getID()
                          << ") has not the good size.");

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

  UInt nb_degree_of_freedom = dof_data.dof->getNbComponent();

  const Array<UInt> & connectivity =
      this->mesh->getConnectivity(type, ghost_type);
  auto conn_begin = connectivity.begin(nb_nodes_per_element);
  auto conn_it = conn_begin;
  auto size_mat = nb_nodes_per_element * nb_degree_of_freedom;

  Vector<Int> element_eq_nb(nb_degree_of_freedom * nb_nodes_per_element);
  auto el_mat_it = elementary_mat.begin(size_mat, size_mat);

  for (UInt e = 0; e < nb_element; ++e, ++el_mat_it) {
    if (filter_it) {
      conn_it = conn_begin + *filter_it;
    }

    this->extractElementEquationNumber(equation_number, *conn_it,
                                       nb_degree_of_freedom, element_eq_nb);
    std::transform(element_eq_nb.begin(), element_eq_nb.end(),
                   element_eq_nb.begin(), [&](auto && local) {
                     return this->localToGlobalEquationNumber(local);
                   });

    if (filter_it) {
      ++filter_it;
    } else {
      ++conn_it;
    }

    A.addValues(element_eq_nb, element_eq_nb, *el_mat_it,
                elemental_matrix_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <typename Mat>
void DOFManager::assemblePreassembledMatrix_(Mat & A, const ID & dof_id_m,
                                             const ID & dof_id_n,
                                             const TermsToAssemble & terms) {
  const auto & equation_number_m = this->getLocalEquationsNumbers(dof_id_m);
  const auto & equation_number_n = this->getLocalEquationsNumbers(dof_id_n);

  for (const auto & term : terms) {
    auto gi = this->localToGlobalEquationNumber(equation_number_m(term.i()));
    auto gj = this->localToGlobalEquationNumber(equation_number_n(term.j()));
    A.add(gi, gj, term);
  }
}
/* -------------------------------------------------------------------------- */

} // namespace akantu

#endif /* AKANTU_DOF_MANAGER_INLINE_IMPL_HH_ */
