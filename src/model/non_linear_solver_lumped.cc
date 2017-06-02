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
#include "static_communicator.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

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

  const Array<Real> & A = this->dof_manager.getLumpedMatrix("M");
  Array<Real> & x = this->dof_manager.getGlobalSolution();
  const Array<Real> & b = this->dof_manager.getResidual();

  x.resize(b.getSize());

  this->dof_manager.updateGlobalBlockedDofs();
  const Array<bool> & blocked_dofs = this->dof_manager.getGlobalBlockedDOFs();

  solver_callback.assembleResidual();

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
  Array<Real>::const_scalar_iterator A_it = A.begin();
  Array<Real>::scalar_iterator x_it = x.begin();
  Array<Real>::scalar_iterator x_end = x.end();
  Array<Real>::const_scalar_iterator b_it = b.begin();

  Array<bool>::const_scalar_iterator blocked_it = blocked_dofs.begin();

  for (; x_it != x_end; ++x_it, ++b_it, ++A_it, ++blocked_it) {
    if (!(*blocked_it)) {
      *x_it = alpha *(*b_it / *A_it);
    }
  }
}

/* -------------------------------------------------------------------------- */

__END_AKANTU__
