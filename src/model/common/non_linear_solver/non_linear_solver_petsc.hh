/**
 * @file   non_linear_solver_petsc.hh
 *
 * @author Nicolas Richart
 *
 * @date creation  Tue Jan 01 2019
 *
 * @brief A Documented file.
 *
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
/* -------------------------------------------------------------------------- */
#include <petscsnes.h>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_NON_LINEAR_SOLVER_PETSC_HH_
#define AKANTU_NON_LINEAR_SOLVER_PETSC_HH_

namespace akantu {
class DOFManagerPETSc;
class NonLinearSolverPETScCallback;
class SolverVectorPETSc;
} // namespace akantu

namespace akantu {

class NonLinearSolverPETSc : public NonLinearSolver {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  NonLinearSolverPETSc(DOFManagerPETSc & dof_manager,
                       const NonLinearSolverType & non_linear_solver_type,
                       const ID & id = "non_linear_solver_petsc");

  ~NonLinearSolverPETSc() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// solve the system described by the jacobian matrix, and rhs contained in
  /// the dof manager
  void solve(SolverCallback & callback) override;

  /// parse the arguments from the input file
  void parseSection(const ParserSection & section) override;
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  static PetscErrorCode FormFunction(SNES snes, Vec dx, Vec f, void * ctx);
  static PetscErrorCode FormJacobian(SNES snes, Vec dx, Mat J, Mat P,
                                     void * ctx);

  void set_param(const ID & param, const std::string & value) override;

  DOFManagerPETSc & dof_manager;

  /// PETSc non linear solver
  SNES snes;
  SNESConvergedReason reason;

  SolverCallback * callback{nullptr};

  std::unique_ptr<SolverVectorPETSc> x;
  std::unique_ptr<NonLinearSolverPETScCallback> ctx;

  Int n_iter{0};
};

namespace debug {
  class SNESNotConvergedException : public NLSNotConvergedException {
  public:
    SNESNotConvergedException(SNESConvergedReason reason, UInt niter,
                              Real error, Real absolute_tolerance,
                              Real relative_tolerance, UInt max_iterations)
        : NLSNotConvergedException(relative_tolerance, niter, error),
          reason(reason), absolute_tolerance(absolute_tolerance),
          max_iterations(max_iterations) {}
    SNESConvergedReason reason;
    Real absolute_tolerance;
    UInt max_iterations;
  };
} // namespace debug

} // namespace akantu

#endif /* AKANTU_NON_LINEAR_SOLVER_PETSC_HH_ */
