/**
 * @file   non_linear_solver_newton_raphson.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Wed Jan 31 2018
 *
 * @brief  Default implementation of NonLinearSolver, in case no external
 * library
 * is there to do the job
 *
 * @section LICENSE
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
#include "non_linear_solver.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_NON_LINEAR_SOLVER_NEWTON_RAPHSON_HH__
#define __AKANTU_NON_LINEAR_SOLVER_NEWTON_RAPHSON_HH__

namespace akantu {
class DOFManagerDefault;
class SparseSolverMumps;
class SolverVector;
} // namespace akantu

namespace akantu {

class NonLinearSolverNewtonRaphson : public NonLinearSolver {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  NonLinearSolverNewtonRaphson(
      DOFManagerDefault & dof_manager,
      const NonLinearSolverType & non_linear_solver_type,
      const ID & id = "non_linear_solver_newton_raphson", UInt memory_id = 0);
  ~NonLinearSolverNewtonRaphson() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// Function that solve the non linear system described by the dof manager and
  /// the solver callback functions
  void solve(SolverCallback & solver_callback) override;

  AKANTU_GET_MACRO_NOT_CONST(Solver, *solver, SparseSolverMumps &);
  AKANTU_GET_MACRO(Solver, *solver, const SparseSolverMumps &);

protected:
  /// test the convergence compare norm of array to convergence_criteria
  bool testConvergence(const SolverVector & array);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  DOFManagerDefault & dof_manager;

  /// Sparse solver used for the linear solves
  std::unique_ptr<SparseSolverMumps> solver;

  /// Type of convergence criteria
  SolveConvergenceCriteria convergence_criteria_type;

  /// convergence threshold
  Real convergence_criteria;

  /// Max number of iterations
  int max_iterations;

  /// Number of iterations at last solve call
  int n_iter{0};

  /// Convergence error at last solve call
  Real error{0.};

  /// Did the last call to solve reached convergence
  bool converged{false};

  /// Force a re-computation of the jacobian matrix
  bool force_linear_recompute{true};
};

} // namespace akantu

#endif /* __AKANTU_NON_LINEAR_SOLVER_NEWTON_RAPHSON_HH__ */
