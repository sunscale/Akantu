/**
 * @file   non_linear_solver_default.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Tue Aug 25 00:57:00 2015
 *
 * @brief  Implementation of the default NonLinearSolver
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

__BEGIN_AKANTU__

NonLinearSolverDefault::NonLinearSolverDefault(
    DOFManagerDefault & dof_manager,
    const NonLinearSolverType & non_linear_solver_type, const ID & id,
    UInt memory_id)
    : NonLinearSolver(dof_manager, non_linear_solver_type, id, memory_id) {}

/* -------------------------------------------------------------------------- */
NonLinearSolverDefault::~NonLinearSolverDefault() {}

/* ------------------------------------------------------------------------ */
void NonLinearSolverDefault::solve() {
  EventManager::sendEvent(NonLinearSolver::BeforeNonLinearSolverSolve(method));

  switch (this->non_linear_solver_type) {
  case _nls_linear:
  case _nls_newton_raphson:
    break;
  case _nls_newton_raphson_modified:
    this->solver_callback.assembleJacobian();
    break;
  default:
    AKANTU_DEBUG_ERROR("The resolution method " << cmethod << " has not been implemented!");
  }

  this->n_iter = 0;
  bool converged = false;
  Real error = 0.;

  if(this->criteria == _scc_residual) {
    converged = this->testConvergence<_scc_residual>(this->convergence_criteria, error);
    if(converged) return converged;
  }

  do {
    if (this->non_linear_solver_type == _nlsnewton_raphson)
      this->solver_callback.assembleJacobian();

    this->solver.solve();

    EventManager::sendEvent(NonLinearSolver::AfterSparseSolve(method));

    if(criteria == _scc_residual) this->solver_callback.assembleResidual();

    converged = this->testConvergence<criteria> (this->convergence_criteria, error);

    if(criteria == _scc_solution && !converged) this->solver_callback.assembleResidual();
    //this->dump();

    this->n_iter++;
    AKANTU_DEBUG_INFO("[" << criteria << "] Convergence iteration "
		      << std::setw(std::log10(max_iteration)) << this->n_iter
		      << ": error " << error << (converged ? " < " : " > ") << tolerance);

  } while (!converged && this->n_iter < this->max_iterations);

  // this makes sure that you have correct strains and stresses after the solveStep function (e.g., for dumping)
  if(criteria == _scc_solution) this->updateResidual();

  if (converged) {
    EventManager::sendEvent(SolidMechanicsModelEvent::AfterNonLinearSolverSolves(method));
  } else if(this->n_iter == max_iteration) {
    AKANTU_DEBUG_WARNING("[" << criteria << "] Convergence not reached after "
			 << std::setw(std::log10(max_iteration)) << this->n_iter <<
			 " iteration" << (this->n_iter == 1 ? "" : "s") << "!" << std::endl);
  }

  return converged;


}

/* -------------------------------------------------------------------------- */
void NonLinearSolverDefault::setParameters(
    const ParserSection & parameters_section) {}

/* -------------------------------------------------------------------------- */

__END_AKANTU__
