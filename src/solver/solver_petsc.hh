/**
 * @file   solver_petsc.hh
 *
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue May 13 2014
 * @date last modification: Mon Jun 19 2017
 *
 * @brief  Solver class interface for the petsc solver
 *
 *
 * Copyright (©) 2014-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "sparse_solver.hh"
/* -------------------------------------------------------------------------- */
#include <petscksp.h>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_SOLVER_PETSC_HH_
#define AKANTU_SOLVER_PETSC_HH_

namespace akantu {
class SparseMatrixPETSc;
class DOFManagerPETSc;
} // namespace akantu

namespace akantu {

class SolverPETSc : public SparseSolver {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  SolverPETSc(DOFManagerPETSc & dof_manager, const ID & matrix_id,
              const ID & id = "solver_petsc", const MemoryID & memory_id = 0);

  ~SolverPETSc() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// create the solver context and set the matrices
  virtual void setOperators();
  void solve() override;

private:
  /// DOFManager correctly typed
  DOFManagerPETSc & dof_manager;

  /// PETSc linear solver
  KSP ksp;

  /// Matrix defining the system of equations
  SparseMatrixPETSc & matrix;
};


} // namespace akantu

#endif /* AKANTU_SOLVER_PETSC_HH_ */
