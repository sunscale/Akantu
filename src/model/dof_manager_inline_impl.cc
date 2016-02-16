/**
 * @file   dof_manager_inline_impl.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Wed Aug 12 11:07:01 2015
 *
 * @brief
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

#ifndef __AKANTU_DOF_MANAGER_INLINE_IMPL_CC__
#define __AKANTU_DOF_MANAGER_INLINE_IMPL_CC__

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
inline void DOFManager::extractElementEquationNumber(
    const Array<UInt> & equation_numbers, const Vector<UInt> & connectivity,
    UInt nb_degree_of_freedom, Vector<UInt> & element_equation_number) {
  for (UInt i = 0, ld = 0; i < connectivity.size(); ++i) {
    UInt n = connectivity(i);
    for (UInt d = 0; d < nb_degree_of_freedom; ++d, ++ld) {
      element_equation_number(ld) =
          equation_numbers(n * nb_degree_of_freedom + d);
    }
  }
}

/* -------------------------------------------------------------------------- */
template <class S>
inline void DOFManager::localToGlobalEquationNumber(S & inout) {
  for (UInt i = 0; i < inout.size(); ++i) {
    inout(i) = this->global_equation_number(inout(i));
  }
}

template<>
inline void DOFManager::localToGlobalEquationNumber<UInt>(UInt & inout) {
  inout = this->global_equation_number(inout);
}

/* -------------------------------------------------------------------------- */
inline UInt DOFManager::globalToLocalEquationNumber(UInt global) const {
  equation_numbers_map::const_iterator it =
      this->global_to_local_mapping.find(global);
  AKANTU_DEBUG_ASSERT(it != this->global_to_local_mapping.end(),
                      "This global equation number "
                          << global << " does not exists in " << this->id);

  return it->second;
}
/* -------------------------------------------------------------------------- */
inline DOFManager::DOFData & DOFManager::getDOFData(const ID & dof_id) {
  DOFStorage::iterator it = this->dofs.find(dof_id);
  if (it == this->dofs.end()) {
    AKANTU_EXCEPTION("The dof " << dof_id << " does not exists in "
                                << this->id);
  }
  return *it->second;
}

/* -------------------------------------------------------------------------- */
const DOFManager::DOFData & DOFManager::getDOFData(const ID & dof_id) const {
  DOFStorage::const_iterator it = this->dofs.find(dof_id);
  if (it == this->dofs.end()) {
    AKANTU_EXCEPTION("The dof " << dof_id << " does not exists in "
                                << this->id);
  }
  return *it->second;
}

/* -------------------------------------------------------------------------- */
inline const Array<UInt> &
DOFManager::getLocalEquationNumbers(const ID & dof_id) const {
  return this->getDOFData(dof_id).local_equation_number;
}

/* -------------------------------------------------------------------------- */
inline Array<Real> & DOFManager::getDOFs(const ID & dofs_id) {
  return *(this->getDOFData(dofs_id).dof);
}

/* -------------------------------------------------------------------------- */
inline Array<Real> & DOFManager::getDOFsDerivatives(const ID & dofs_id,
                                                    UInt order) {
  std::vector<Array<Real> *> & derivatives =
      this->getDOFData(dofs_id).dof_derivatives;
  if (order >= derivatives.size())
    AKANTU_EXCEPTION("No derivatives of order " << order << " present in "
                                                << this->id << " for dof "
                                                << dofs_id);

  return *derivatives[order];
}

/* -------------------------------------------------------------------------- */
inline const Array<Real> & DOFManager::getSolution(const ID & dofs_id) const {
  return *(this->getDOFData(dofs_id).solution);
}

/* -------------------------------------------------------------------------- */
inline const Array<bool> &
DOFManager::getBlockedDOFs(const ID & dofs_id) const {
  return *(this->getDOFData(dofs_id).blocked_dofs);
}

/* -------------------------------------------------------------------------- */
__END_AKANTU__

#endif /* __AKANTU_DOF_MANAGER_INLINE_IMPL_CC__ */
