/**
 * @file   node_synchronizer_inline_impl.hh
 *
 * @author Nicolas Richart
 *
 * @date creation  mar jan 14 2020
 *
 * @brief A Documented file.
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
#include "node_synchronizer.hh"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_NODE_SYNCHRONIZER_INLINE_IMPL_HH__
#define __AKANTU_NODE_SYNCHRONIZER_INLINE_IMPL_HH__

namespace akantu {

/* -------------------------------------------------------------------------- */
inline UInt NodeSynchronizer::canScatterSize() {
  return mesh.getNbNodes();
}

/* -------------------------------------------------------------------------- */
inline UInt NodeSynchronizer::gatheredSize() {
  return mesh.getNbGlobalNodes();
}

/* -------------------------------------------------------------------------- */
inline UInt NodeSynchronizer::localToGlobalEntity(const UInt & local) {
  return mesh.getNodeGlobalId(local);
}

} // akantu

#endif // __AKANTU_NODE_SYNCHRONIZER_INLINE_IMPL_HH__
