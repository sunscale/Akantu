/**
 * @file   non_linear_solver_lumped.hh
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
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_NON_LINEAR_SOLVER_LUMPED_HH__
#define __AKANTU_NON_LINEAR_SOLVER_LUMPED_HH__

namespace akantu {
class DOFManagerDefault;
}

namespace akantu {

class NonLinearSolverLumped : public NonLinearSolver {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  NonLinearSolverLumped(DOFManagerDefault & dof_manager,
                        const NonLinearSolverType & non_linear_solver_type,
                        const ID & id = "non_linear_solver_lumped",
                        UInt memory_id = 0);
  ~NonLinearSolverLumped() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// Function that solve the non linear system described by the dof manager and
  /// the solver callback functions
  void solve(SolverCallback & solver_callback) override;

  static void solveLumped(const Array<Real> & A, Array<Real> & x,
                          const Array<Real> & b, Real alpha,
                          const Array<bool> & blocked_dofs);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  DOFManagerDefault & dof_manager;

  /// Coefficient to apply between x and A^{-1} b
  Real alpha;
};

} // namespace akantu

#endif /* __AKANTU_NON_LINEAR_SOLVER_LUMPED_HH__ */
