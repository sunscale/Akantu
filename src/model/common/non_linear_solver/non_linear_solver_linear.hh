/**
 * @file   non_linear_solver_linear.hh
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
#include "sparse_solver_mumps.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_NON_LINEAR_SOLVER_LINEAR_HH_
#define AKANTU_NON_LINEAR_SOLVER_LINEAR_HH_

namespace akantu {
class DOFManagerDefault;
}

namespace akantu {

class NonLinearSolverLinear : public NonLinearSolver {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  NonLinearSolverLinear(DOFManagerDefault & dof_manager,
                        const NonLinearSolverType & non_linear_solver_type,
                        const ID & id = "non_linear_solver_linear");
  ~NonLinearSolverLinear() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// Function that solve the non linear system described by the dof manager and
  /// the solver callback functions
  void solve(SolverCallback & solver_callback) override;

  AKANTU_GET_MACRO_NOT_CONST(Solver, solver, SparseSolverMumps &);
  AKANTU_GET_MACRO(Solver, solver, const SparseSolverMumps &);
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  DOFManagerDefault & dof_manager;

  /// Sparse solver used for the linear solves
  SparseSolverMumps solver;
};

} // namespace akantu

#endif /* AKANTU_NON_LINEAR_SOLVER_LINEAR_HH_ */
