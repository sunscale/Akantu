/**
 * @file   dof_manager_inline_impl.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Feb 21 2013
 * @date last modification: Wed Jan 31 2018
 *
 * @brief  inline functions of the dof manager
 *
 * @section LICENSE
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
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_DOF_MANAGER_INLINE_IMPL_CC__
#define __AKANTU_DOF_MANAGER_INLINE_IMPL_CC__

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
const Array<UInt> & DOFManager::getEquationsNumbers(const ID & dof_id) const {
  return getDOFData(dof_id).local_equation_number;
}

/* -------------------------------------------------------------------------- */
template <class _DOFData>
inline _DOFData & DOFManager::getDOFDataTyped(const ID & dof_id) {
  return dynamic_cast<_DOFData &>(this->getDOFData(dof_id));
}

/* -------------------------------------------------------------------------- */
template <class _DOFData>
inline const _DOFData & DOFManager::getDOFDataTyped(const ID & dof_id) const {
  return dynamic_cast<const _DOFData &>(this->getDOFData(dof_id));
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
  std::vector<Array<Real> *> & derivatives =
      this->getDOFData(dofs_id).dof_derivatives;
  if ((order > derivatives.size()) || (derivatives[order - 1] == nullptr))
    AKANTU_EXCEPTION("No derivatives of order " << order << " present in "
                                                << this->id << " for dof "
                                                << dofs_id);

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
inline const Array<bool> &
DOFManager::getBlockedDOFs(const ID & dofs_id) const {
  return *(this->getDOFData(dofs_id).blocked_dofs);
}

/* -------------------------------------------------------------------------- */
inline bool DOFManager::hasBlockedDOFs(const ID & dofs_id) const {
  return (this->getDOFData(dofs_id).blocked_dofs != nullptr);
}

/* -------------------------------------------------------------------------- */

} // akantu

#endif /* __AKANTU_DOF_MANAGER_INLINE_IMPL_CC__ */
