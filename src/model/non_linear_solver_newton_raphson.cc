/**
 * @file   non_linear_solver_newton_raphson.cc
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
#include "non_linear_solver_newton_raphson.hh"
#include "dof_manager_default.hh"
#include "solver_callback.hh"
#include "static_communicator.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
NonLinearSolverNewtonRaphson::NonLinearSolverNewtonRaphson(
    DOFManagerDefault & dof_manager,
    const NonLinearSolverType & non_linear_solver_type, const ID & id,
    UInt memory_id)
    : NonLinearSolver(dof_manager, non_linear_solver_type, id, memory_id),
      dof_manager(dof_manager),
      solver(dof_manager, "J", id + ":sparse_solver", memory_id), n_iter(0),
      error(0.), converged(false) {

  this->supported_type.insert(_nls_newton_raphson_modified);
  this->supported_type.insert(_nls_newton_raphson);
  this->supported_type.insert(_nls_linear);

  this->checkIfTypeIsSupported();

  this->registerParam("threshold", convergence_criteria, 1e-10, _pat_parsmod,
                      "Threshold to consider results as converged");
  this->registerParam("convergence_type", convergence_criteria_type,
                      _scc_solution, _pat_parsmod,
                      "Type of convergence criteria");
  this->registerParam("max_iterations", max_iterations, UInt(10), _pat_parsmod,
                      "Max number of iterations");
  this->registerParam("error", error, _pat_readable, "Last reached error");
  this->registerParam("nb_iterations", n_iter, _pat_readable,
                      "Last reached number of iterations");
  this->registerParam("converged", converged, _pat_readable,
                      "Did last solve converged");
  this->registerParam("force_linear_recompute", force_linear_recompute, true,
                      _pat_modifiable,
                      "Force reassembly of the jacobian matrix");
}

/* -------------------------------------------------------------------------- */
NonLinearSolverNewtonRaphson::~NonLinearSolverNewtonRaphson() {}

/* ------------------------------------------------------------------------ */
void NonLinearSolverNewtonRaphson::solve(SolverCallback & solver_callback) {
  solver_callback.predictor();

  solver_callback.assembleResidual();

  if (this->non_linear_solver_type == _nls_newton_raphson_modified ||
      (this->non_linear_solver_type == _nls_linear &&
       this->force_linear_recompute)) {
    solver_callback.assembleJacobian();
    this->force_linear_recompute = false;
  }

  this->n_iter = 0;
  this->converged = false;

  if (this->convergence_criteria_type == _scc_residual) {
    this->converged = this->testConvergence(this->dof_manager.getResidual());

    if (this->converged)
      return;
  }

  do {
    if (this->non_linear_solver_type == _nls_newton_raphson)
      solver_callback.assembleJacobian();

    this->solver.solve();

    solver_callback.corrector();

    // EventManager::sendEvent(NonLinearSolver::AfterSparseSolve(method));

    if (this->convergence_criteria_type == _scc_residual) {
      solver_callback.assembleResidual();
      this->converged = this->testConvergence(this->dof_manager.getResidual());
    } else {
      this->converged =
          this->testConvergence(this->dof_manager.getGlobalSolution());
    }

    if (this->convergence_criteria_type == _scc_solution && !this->converged)
      solver_callback.assembleResidual();
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
    solver_callback.assembleResidual();

  if (this->converged) {
    //    EventManager::sendEvent(
    //   SolidMechanicsModelEvent::AfterNonLinearSolverSolves(method));
  } else if (this->n_iter == this->max_iterations) {
    AKANTU_CUSTOM_EXCEPTION(debug::NLSNotConvergedException(
        this->convergence_criteria, this->n_iter));

    AKANTU_DEBUG_WARNING("[" << this->convergence_criteria_type
                             << "] Convergence not reached after "
                             << std::setw(std::log10(this->max_iterations))
                             << this->n_iter << " iteration"
                             << (this->n_iter == 1 ? "" : "s") << "!");
  }

  return;
}

/* -------------------------------------------------------------------------- */
bool NonLinearSolverNewtonRaphson::testConvergence(const Array<Real> & array) {
  AKANTU_DEBUG_IN();

  const Array<bool> & blocked_dofs = this->dof_manager.getGlobalBlockedDOFs();

  UInt nb_degree_of_freedoms = array.getSize();

  Array<Real>::const_scalar_iterator arr_it = array.begin();
  Array<bool>::const_scalar_iterator bld_it = blocked_dofs.begin();

  Real norm = 0.;
  for (UInt n = 0; n < nb_degree_of_freedoms; ++n, ++arr_it, ++bld_it) {
    bool is_local_node = this->dof_manager.isLocalOrMasterDOF(n);
    if ((*bld_it == 0) && is_local_node) {
      norm += *arr_it * *arr_it;
    }
  }

  StaticCommunicator::getStaticCommunicator().allReduce(norm, _so_sum);

  norm = std::sqrt(norm);

  AKANTU_DEBUG_ASSERT(!Math::isnan(norm),
                      "Something went wrong in the solve phase");

  this->error = norm;

  return (error < this->convergence_criteria);
}

/* -------------------------------------------------------------------------- */

__END_AKANTU__
