/**
 * @file   non_linear_solver_petsc.cc
 *
 * @author Nicolas Richart
 *
 * @date creation  Mon Dec 31 2018
 *
 * @brief A Documented file.
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
#include "non_linear_solver_petsc.hh"
#include "dof_manager_petsc.hh"
#include "mpi_communicator_data.hh"
/* -------------------------------------------------------------------------- */
#include <petscoptions.h>
/* -------------------------------------------------------------------------- */

namespace akantu {

NonLinearSolverPETSc::NonLinearSolverPETSc(
    DOFManagerPETSc & dof_manager,
    const NonLinearSolverType & non_linear_solver_type, const ID & id,
    UInt memory_id)
    : NonLinearSolver(dof_manager, non_linear_solver_type, id, memory_id),
      dof_manager(dof_manager) {
  std::unordered_map<NonLinearSolverType, SNESType>
      petsc_non_linear_solver_types{
          {NonLinearSolverType::_newton_raphson, SNESNEWTONLS},
          {NonLinearSolverType::_linear, SNESKSPONLY},
          {NonLinearSolverType::_gmres, SNESNGMRES},
          {NonLinearSolverType::_bfgs, SNESQN},
          {NonLinearSolverType::_cg, SNESNCG}};

  for (const auto & pair : petsc_non_linear_solver_types) {
    supported_type.insert(pair.first);
  }

  this->checkIfTypeIsSupported();

  auto mpi_comm = dof_manager.getMPIComm();

  PETSc_call(SNESCreate, mpi_comm, &snes);

  auto it = petsc_non_linear_solver_types.find(non_linear_solver_type);
  if (it != petsc_non_linear_solver_types.end()) {
    PETSc_call(SNESSetType, snes, it->second);
  }
}

/* -------------------------------------------------------------------------- */
NonLinearSolverPETSc::~NonLinearSolverPETSc() {
  PETSc_call(SNESDestroy, &snes);
}

/* -------------------------------------------------------------------------- */
void NonLinearSolverPETSc::solve(SolverCallback & /*callback*/) {

}

/* -------------------------------------------------------------------------- */
void NonLinearSolverPETSc::parseSection(const ParserSection & section) {
  auto parameters = section.getParameters();
  for (auto && param : range(parameters.first, parameters.second)) {
    PetscOptionsSetValue(NULL, param.getName().c_str(), param.getValue().c_str());
  }

  SNESSetFromOptions(snes);

  PetscOptionsClear(NULL);
}

} // namespace akantu
