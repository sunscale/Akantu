/**
 * @file   non_linear_solver_lumped.cc
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
#include "non_linear_solver_lumped.hh"
#include "dof_manager_default.hh"
#include "solver_callback.hh"
#include "communicator.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
NonLinearSolverLumped::NonLinearSolverLumped(
    DOFManagerDefault & dof_manager,
    const NonLinearSolverType & non_linear_solver_type, const ID & id,
    UInt memory_id)
    : NonLinearSolver(dof_manager, non_linear_solver_type, id, memory_id),
      dof_manager(dof_manager) {
  this->supported_type.insert(_nls_lumped);
  this->checkIfTypeIsSupported();

  this->registerParam("b_a2x", this->alpha, 1., _pat_parsmod,
                      "Conversion coefficient between x and A^{-1} b");
}

/* -------------------------------------------------------------------------- */
NonLinearSolverLumped::~NonLinearSolverLumped() {}

/* ------------------------------------------------------------------------ */
void NonLinearSolverLumped::solve(SolverCallback & solver_callback) {
  this->dof_manager.updateGlobalBlockedDofs();
  solver_callback.predictor();

  auto & x = this->dof_manager.getGlobalSolution();
  const auto & b = this->dof_manager.getResidual();

  x.resize(b.size());

  this->dof_manager.updateGlobalBlockedDofs();
  const Array<bool> & blocked_dofs = this->dof_manager.getGlobalBlockedDOFs();

  solver_callback.assembleResidual();

  const auto & A = this->dof_manager.getLumpedMatrix("M");
  // alpha is the conversion factor from from force/mass to acceleration needed
  // in model coupled with atomistic \todo find a way to define alpha per dof
  // type

  this->solveLumped(A, x, b, blocked_dofs, alpha);
  this->dof_manager.splitSolutionPerDOFs();

  solver_callback.corrector();
}

/* -------------------------------------------------------------------------- */
void NonLinearSolverLumped::solveLumped(const Array<Real> & A, Array<Real> & x,
                                        const Array<Real> & b,
                                        const Array<bool> & blocked_dofs,
                                        Real alpha) {
  auto A_it = A.begin();
  auto x_it = x.begin();
  auto x_end = x.end();
  auto b_it = b.begin();
  auto blocked_it = blocked_dofs.begin();

  for (; x_it != x_end; ++x_it, ++b_it, ++A_it, ++blocked_it) {
    if (!(*blocked_it)) {
      *x_it = alpha *(*b_it / *A_it);
    }
  }
}

/* -------------------------------------------------------------------------- */

} // akantu
