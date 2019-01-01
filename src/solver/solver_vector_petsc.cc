/**
 * @file   solver_vector_petsc.cc
 *
 * @author Nicolas Richart
 *
 * @date creation  Tue Jan 01 2019
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
#include "solver_vector_petsc.hh"
#include "dof_manager_petsc.hh"
#include "mpi_communicator_data.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
SolverVectorPETSc::SolverVectorPETSc(DOFManagerPETSc & dof_manager,
                                     const ID & id)
    : SolverVector(dof_manager, id), dof_manager(dof_manager) {
  auto && mpi_comm = dof_manager.getMPIComm();
  PETSc_call(VecCreate, mpi_comm, &vector);
  PETSc_call(PetscObjectSetName, reinterpret_cast<PetscObject>(vector), id.c_str());
}

/* -------------------------------------------------------------------------- */
SolverVectorPETSc::SolverVectorPETSc(const SolverVectorPETSc & vector,
                                     const ID & id)
    : SolverVector(vector, id), dof_manager(vector.dof_manager) {
  PETSc_call(VecDuplicate, vector.vector, &this->vector);
  PETSc_call(PetscObjectSetName, reinterpret_cast<PetscObject>(this->vector), id.c_str());

  PETSc_call(VecCopy, vector.vector, this->vector);
}

/* -------------------------------------------------------------------------- */
SolverVectorPETSc::~SolverVectorPETSc() {
  if (vector) {
    PETSc_call(VecDestroy, &vector);
  }
}

/* -------------------------------------------------------------------------- */
void SolverVectorPETSc::resize() {
  PETSc_call(VecSetSizes, vector, dof_manager.getPureLocalSystemSize(),
             dof_manager.getSystemSize());

  auto & is_ltog_mapping = dof_manager.getISLocalToGlobalMapping();

  PETSc_call(VecSetLocalToGlobalMapping, vector, is_ltog_mapping);
}

/* -------------------------------------------------------------------------- */
void SolverVectorPETSc::clear() {
  PETSc_call(VecSet, vector, 0.);
  applyModifications();
}

/* -------------------------------------------------------------------------- */
void SolverVectorPETSc::applyModifications() {
  PETSc_call(VecAssemblyBegin, vector);
  PETSc_call(VecAssemblyEnd, vector);
}

} // namespace akantu
