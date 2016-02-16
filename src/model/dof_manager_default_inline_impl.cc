/**
 * @file   dof_manager_default_inline_impl.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Wed Aug 12 10:57:47 2015
 *
 * @brief  Implementation of the DOFManagerDefault inline functions
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

#ifndef __AKANTU_DOF_MANAGER_DEFAULT_INLINE_IMPL_CC__
#define __AKANTU_DOF_MANAGER_DEFAULT_INLINE_IMPL_CC__

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
inline bool DOFManagerDefault::isLocalOrMasterDOF(UInt dof_num) {
  Int dof_type = this->dofs_type(dof_num);
  return (dof_type == -2) || (dof_type == -1);
}

/* -------------------------------------------------------------------------- */

__END_AKANTU__

#endif /* __AKANTU_DOF_MANAGER_DEFAULT_INLINE_IMPL_CC__ */
