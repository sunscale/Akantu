/**
 * @file   non_local_manager_callback.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jul 21 2017
 * @date last modification: Tue Sep 19 2017
 *
 * @brief  Callback functions for the non local manager
 *
 *
 * Copyright (©) 2016-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_common.hh"
#include "element_type_map.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_NON_LOCAL_MANAGER_CALLBACK_HH_
#define AKANTU_NON_LOCAL_MANAGER_CALLBACK_HH_

namespace akantu {
class NonLocalManager;
} // namespace akantu

namespace akantu {

class NonLocalManagerCallback {
public:
  virtual void initializeNonLocal() {}

  /* ------------------------------------------------------------------------ */
  virtual void
  insertIntegrationPointsInNeighborhoods(GhostType ghost_type) = 0;

  virtual void computeNonLocalStresses(GhostType ghost_type) = 0;

  /// update the values of the non local internal
  virtual void updateLocalInternal(ElementTypeMapReal & internal_flat,
                                   GhostType ghost_type,
                                   ElementKind kind) = 0;

  /// copy the results of the averaging in the materials
  virtual void updateNonLocalInternal(ElementTypeMapReal & internal_flat,
                                      GhostType ghost_type,
                                      ElementKind kind) = 0;
};

} // namespace akantu

#endif /* AKANTU_NON_LOCAL_MANAGER_CALLBACK_HH_ */
