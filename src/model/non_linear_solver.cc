/**
 * @file   non_linear_solver.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Tue Oct 13 15:34:43 2015
 *
 * @brief  Implementation of the base class NonLinearSolver
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
#include "non_linear_solver.hh"
#include "solver_callback.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
NonLinearSolver::NonLinearSolver(
    DOFManager & dof_manager,
    const NonLinearSolverType & non_linear_solver_type, const ID & id,
    UInt memory_id)
    : Memory(id, memory_id), Parsable(_st_non_linear_solver, id),
      _dof_manager(dof_manager),
      non_linear_solver_type(non_linear_solver_type) {

  this->registerParam("type", this->non_linear_solver_type, _pat_parsable,
                      "Non linear solver type");
}

/* -------------------------------------------------------------------------- */
NonLinearSolver::~NonLinearSolver() {}

/* -------------------------------------------------------------------------- */
void NonLinearSolver::checkIfTypeIsSupported() {
  if (this->supported_type.find(this->non_linear_solver_type) ==
      this->supported_type.end()) {
    AKANTU_EXCEPTION("The resolution method "
                     << this->non_linear_solver_type
                     << " is not implemented in the non linear solver "
                     << this->id << "!");
  }
}

/* -------------------------------------------------------------------------- */

} // akantu
