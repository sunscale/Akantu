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
#include "non_linear_solver_default.hh"
#include "dof_manager_default.hh"
#include "non_linear_solver_callback.hh"
#include "static_communicator.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
NonLinearSolverDefault::NonLinearSolverDefault(
    DOFManagerDefault & dof_manager,
    const NonLinearSolverType & non_linear_solver_type, const ID & id,
    UInt memory_id)
    : NonLinearSolver(dof_manager, non_linear_solver_type, id, memory_id),
      dof_manager(dof_manager),
      solver(dof_manager, "jacobian", id + ":sparse_solver", memory_id),
      convergence_criteria_type(_scc_solution), convergence_criteria(1e-10),
      max_iterations(10), n_iter(0), error(0.), converged(false) {}

/* -------------------------------------------------------------------------- */
NonLinearSolverDefault::~NonLinearSolverDefault() {}

/* ------------------------------------------------------------------------ */
void NonLinearSolverDefault::solve() {
  // EventManager::sendEvent(NonLinearSolver::BeforeNonLinearSolverSolve(method));
  this->callbackPredictors();

  switch (this->non_linear_solver_type) {
  case _nls_linear:
  case _nls_newton_raphson:
    break;
  case _nls_newton_raphson_modified:
    this->callbackAssembleJacobian();
    break;
  default:
    AKANTU_DEBUG_ERROR("The resolution method "
                       << this->non_linear_solver_type
                       << " has not been implemented!");
  }

  this->n_iter = 0;
  this->converged = false;

  if (this->convergence_criteria_type == _scc_residual) {
    this->converged = this->testConvergence(this->dof_manager.getResidual());

    if (this->converged)  return;
  }

  do {
    if (this->non_linear_solver_type == _nls_newton_raphson)
      this->callbackAssembleJacobian();

    this->solver.solve();

    this->callbackCorrectors();

    // EventManager::sendEvent(NonLinearSolver::AfterSparseSolve(method));

    if (this->convergence_criteria_type == _scc_residual) {
      this->callbackAssembleResidual();
      this->converged = this->testConvergence(this->dof_manager.getResidual());
    } else {
      this->converged = this->testConvergence(this->dof_manager.getGlobalSolution());
    }

    if (this->convergence_criteria_type == _scc_solution && !this->converged)
      this->callbackAssembleResidual();
    // this->dump();

    this->n_iter++;
    AKANTU_DEBUG_INFO(
        "[" << this->convergence_criteria_type << "] Convergence iteration "
            << std::setw(std::log10(this->max_iterations)) << this->n_iter
            << ": error " << this->error << (this->converged ? " < " : " > ")
            << this->convergence_criteria);

  } while (!this->converged && this->n_iter < this->max_iterations);

  // this makes sure that you have correct strains and stresses after the
  // solveStep function (e.g., for dumping)
  if (this->convergence_criteria_type == _scc_solution)
    this->callbackAssembleResidual();

  if (this->converged) {
    //    EventManager::sendEvent(
    //   SolidMechanicsModelEvent::AfterNonLinearSolverSolves(method));
  } else if (this->n_iter == this->max_iterations) {
    AKANTU_DEBUG_WARNING("[" << this->convergence_criteria_type
                             << "] Convergence not reached after "
                             << std::setw(std::log10(this->max_iterations))
                             << this->n_iter << " iteration"
                             << (this->n_iter == 1 ? "" : "s") << "!"
                             << std::endl);
  }

  return;
}

/* -------------------------------------------------------------------------- */
void NonLinearSolverDefault::setParameters(
    const ParserSection & parameters_section) {}

/* -------------------------------------------------------------------------- */
bool NonLinearSolverDefault::testConvergence(const Array<Real> & array) {
  AKANTU_DEBUG_IN();

  const Array<bool> & blocked_dofs = this->dof_manager.getGlobalBlockedDOFs();

  UInt nb_degree_of_freedoms = array.getSize();

  Array<Real>::const_scalar_iterator arr_it = array.begin();
  Array<bool>::const_scalar_iterator bld_it = blocked_dofs.begin();

  Real norm = 0.;
  for (UInt n = 0; n < nb_degree_of_freedoms; ++n, ++arr_it, ++bld_it) {
    bool is_local_node = this->dof_manager.isLocalOrMasterDOF(n);
    if(!(*bld_it) && is_local_node) {
        norm += *arr_it * *arr_it;
    }
  }

  StaticCommunicator::getStaticCommunicator().allReduce(&norm, 1, _so_sum);

  norm = std::sqrt(norm);

  AKANTU_DEBUG_ASSERT(!Math::isnan(norm), "Something goes wrong in the solve phase");

  this->error = norm;

  return (error < this->convergence_criteria);
}

/* -------------------------------------------------------------------------- */



__END_AKANTU__
