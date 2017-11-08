/**
 * @file   solver_callback.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Wed Mar  2 12:47:17 2016
 *
 * @brief  Default behavior of solver_callbacks
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
#include "solver_callback.hh"
#include "dof_manager.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
SolverCallback::SolverCallback(DOFManager & dof_manager)
    : sc_dof_manager(&dof_manager) {}

/* -------------------------------------------------------------------------- */
SolverCallback::SolverCallback() {}

/* -------------------------------------------------------------------------- */
SolverCallback::~SolverCallback() = default;

/* -------------------------------------------------------------------------- */
void SolverCallback::setDOFManager(DOFManager & dof_manager) {
  this->sc_dof_manager = &dof_manager;
}

/* -------------------------------------------------------------------------- */

} // namespace akantu
