/**
 * @file   time_step_solver.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Aug 18 2015
 * @date last modification: Wed Feb 21 2018
 *
 * @brief  Implementation of common part of TimeStepSolvers
 *
 * @section LICENSE
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "time_step_solver.hh"
#include "dof_manager.hh"
#include "non_linear_solver.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
TimeStepSolver::TimeStepSolver(DOFManager & dof_manager,
                               const TimeStepSolverType & type,
                               NonLinearSolver & non_linear_solver,
                               SolverCallback & solver_callback, const ID & id,
                               UInt memory_id)
    : Memory(id, memory_id), SolverCallback(dof_manager),
      _dof_manager(dof_manager), type(type), time_step(0.),
      solver_callback(&solver_callback), non_linear_solver(non_linear_solver) {
  this->registerSubRegistry("non_linear_solver", non_linear_solver);
}

/* -------------------------------------------------------------------------- */
TimeStepSolver::~TimeStepSolver() = default;

/* -------------------------------------------------------------------------- */
void TimeStepSolver::setIntegrationScheme(
    const ID & dof_id, const IntegrationSchemeType & type,
    IntegrationScheme::SolutionType solution_type) {
  this->setIntegrationSchemeInternal(dof_id, type, solution_type);

  for (auto & pair : needed_matrices) {
    auto & mat_type = pair.second;
    const auto & name = pair.first;
    if (mat_type == _mt_not_defined) {
      mat_type = this->solver_callback->getMatrixType(name);
    }

    if (mat_type == _mt_not_defined)
      continue;

    if (not _dof_manager.hasMatrix(name))
      _dof_manager.getNewMatrix(name, mat_type);
  }
}

/* -------------------------------------------------------------------------- */
MatrixType TimeStepSolver::getCommonMatrixType() {
  MatrixType common_type = _mt_not_defined;
  for (auto & pair : needed_matrices) {
    auto & type = pair.second;
    common_type = std::min(common_type, type);
  }

  AKANTU_DEBUG_ASSERT(common_type != _mt_not_defined,
                      "No type defined for the matrices");

  return common_type;
}

/* -------------------------------------------------------------------------- */
void TimeStepSolver::predictor() {
  AKANTU_DEBUG_ASSERT(
      this->solver_callback != nullptr,
      "This function cannot be called if the solver_callback is not set");

  this->solver_callback->predictor();
}

/* -------------------------------------------------------------------------- */
void TimeStepSolver::corrector() {
  AKANTU_DEBUG_ASSERT(
      this->solver_callback != nullptr,
      "This function cannot be called if the solver_callback is not set");

  this->solver_callback->corrector();
}

/* -------------------------------------------------------------------------- */
void TimeStepSolver::beforeSolveStep() {
  AKANTU_DEBUG_ASSERT(
      this->solver_callback != nullptr,
      "This function cannot be called if the solver_callback is not set");

  this->solver_callback->beforeSolveStep();
}

/* -------------------------------------------------------------------------- */
void TimeStepSolver::afterSolveStep(bool converged) {
  AKANTU_DEBUG_ASSERT(
      this->solver_callback != nullptr,
      "This function cannot be called if the solver_callback is not set");

  this->solver_callback->afterSolveStep(converged);
}

/* -------------------------------------------------------------------------- */
void TimeStepSolver::assembleLumpedMatrix(const ID & matrix_id) {
  AKANTU_DEBUG_ASSERT(
      this->solver_callback != nullptr,
      "This function cannot be called if the solver_callback is not set");

  if (not _dof_manager.hasLumpedMatrix(matrix_id))
    _dof_manager.getNewLumpedMatrix(matrix_id);

  this->solver_callback->assembleLumpedMatrix(matrix_id);
}

/* -------------------------------------------------------------------------- */
void TimeStepSolver::assembleMatrix(const ID & matrix_id) {
  AKANTU_DEBUG_ASSERT(
      this->solver_callback != nullptr,
      "This function cannot be called if the solver_callback is not set");

  auto common_type = this->getCommonMatrixType();

  if (matrix_id != "J") {
    auto type = needed_matrices[matrix_id];
    if (type == _mt_not_defined)
      return;

    if (not _dof_manager.hasMatrix(matrix_id)) {
      _dof_manager.getNewMatrix(matrix_id, type);
    }

    this->solver_callback->assembleMatrix(matrix_id);
    return;
  }

  if (not _dof_manager.hasMatrix("J"))
    _dof_manager.getNewMatrix("J", common_type);

  MatrixType type;
  ID name;
  for (auto & pair : needed_matrices) {
    std::tie(name, type) = pair;
    if (type == _mt_not_defined)
      continue;

    this->solver_callback->assembleMatrix(name);
  }
}

/* -------------------------------------------------------------------------- */
void TimeStepSolver::assembleResidual() {
  AKANTU_DEBUG_ASSERT(
      this->solver_callback != nullptr,
      "This function cannot be called if the solver_callback is not set");

  this->_dof_manager.clearResidual();
  this->solver_callback->assembleResidual();
}

/* -------------------------------------------------------------------------- */
void TimeStepSolver::assembleResidual(const ID & residual_part) {
  AKANTU_DEBUG_ASSERT(
      this->solver_callback != nullptr,
      "This function cannot be called if the solver_callback is not set");

  this->solver_callback->assembleResidual(residual_part);
}

/* -------------------------------------------------------------------------- */

} // namespace akantu
