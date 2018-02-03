/* -------------------------------------------------------------------------- */
#include "non_linear_solver.hh"
#include "petscsnes.h"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_NON_LINEAR_SOLVER_PETSC_HH__
#define __AKANTU_NON_LINEAR_SOLVER_PETSC_HH__

namespace akantu {

class NonLinearSolverPETSc : NonLinearSolver {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  NonLinearSolverPETSc(DOFManagerPETSc & dof_manager,
                       const NonLinearSolverType & non_linear_solver_type,
                       const ID & id = "non_linear_solver_petsc",
                       UInt memory_id = 0);

  ~NonLinearSolverPETSc() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// solve the system described by the jacobian matrix, and rhs contained in
  /// the dof manager
  void solve(SolverCallback & callback) override;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  DOFManagerPETSc & dof_manager;

  SNES snes;
};

} // akantu

#endif /* __AKANTU_NON_LINEAR_SOLVER_PETSC_HH__ */
