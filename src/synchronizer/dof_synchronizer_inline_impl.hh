/**
 * @file   dof_synchronizer_inline_impl.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 17 2011
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  DOFSynchronizer inline implementation
 *
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "communication_buffer.hh"
#include "data_accessor.hh"
#include "dof_manager_default.hh"
#include "dof_synchronizer.hh"
/* -------------------------------------------------------------------------- */
#include <map>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_DOF_SYNCHRONIZER_INLINE_IMPL_HH__
#define __AKANTU_DOF_SYNCHRONIZER_INLINE_IMPL_HH__

namespace akantu {

/* -------------------------------------------------------------------------- */
inline UInt DOFSynchronizer::canScatterSize() {
  return dof_manager.getLocalSystemSize();
}

/* -------------------------------------------------------------------------- */
inline UInt DOFSynchronizer::gatheredSize() {
  return dof_manager.getSystemSize();
}

inline UInt DOFSynchronizer::localToGlobalEntity(const UInt & local) {
  return dof_manager.localToGlobalEquationNumber(local);
}
} // namespace akantu

#endif /* __AKANTU_DOF_SYNCHRONIZER_INLINE_IMPL_HH__ */
