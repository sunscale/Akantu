/**
 * @file   boundary_functions.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Fri Jan 04 2013
 * @date last modification: Fri Feb 23 2018
 *
 * @brief  functions for boundaries
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
#include "aka_common.hh"
/* -------------------------------------------------------------------------- */
#include <vector>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_BOUNDARY_FUNCTIONS_HH_
#define AKANTU_BOUNDARY_FUNCTIONS_HH_

namespace akantu {
class SolidMechanicsModel;
}

namespace akantu {

Real integrateResidual(const std::string & sub_boundary_name,
                       const SolidMechanicsModel & model, UInt dir);

/// this is a fix so that all subboundaries exist on all procs
void boundaryFix(Mesh & mesh,
                 const std::vector<std::string> & sub_boundary_names);

} // namespace akantu

#endif /* AKANTU_BOUNDARY_FUNCTIONS_HH_ */
