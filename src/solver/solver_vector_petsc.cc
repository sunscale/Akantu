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
#include <petscvec.h>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
SolverVectorPETSc::SolverVectorPETSc(DOFManagerPETSc & dof_manager,
                                     const ID & id)
    : SolverVector(dof_manager, id), dof_manager(dof_manager) {}

/* -------------------------------------------------------------------------- */
SolverVectorPETSc::SolverVectorPETSc(const SolverVectorPETSc & vector,
                                     const ID & id)
    : SolverVector(vector, id), dof_manager(vector.dof_manager) {
  if (vector.x) {
    PETSc_call(VecDuplicate, vector.x, &x);
    PETSc_call(VecCopy, vector.x, x);
    detail::PETScSetName(x, id);
  }
}

/* -------------------------------------------------------------------------- */
SolverVectorPETSc::SolverVectorPETSc(Vec x, DOFManagerPETSc & dof_manager,
                                     const ID & id)
    : SolverVector(dof_manager, id), dof_manager(dof_manager) {
  PETSc_call(VecDuplicate, x, &this->x);
  PETSc_call(VecCopy, x, this->x);
  detail::PETScSetName(x, id);
}

/* -------------------------------------------------------------------------- */
SolverVectorPETSc::~SolverVectorPETSc() {
  if (x) {
    PETSc_call(VecDestroy, &x);
  }
}

/* -------------------------------------------------------------------------- */
void SolverVectorPETSc::resize() {
  if (not x) {
    auto && mpi_comm = dof_manager.getMPIComm();
    PETSc_call(VecCreate, mpi_comm, &x);
    detail::PETScSetName(x, id);
  }

  auto local_system_size = dof_manager.getLocalSystemSize();
  auto nb_local_dofs = dof_manager.getPureLocalSystemSize();
  PETSc_call(VecSetSizes, x, nb_local_dofs, PETSC_DECIDE);

  VecType vec_type;
  PETSc_call(VecGetType, x, &vec_type);
  if (vec_type == VECMPI) {
    std::vector<PetscInt> idx;
    idx.reserve(local_system_size - nb_local_dofs);

    for (auto && d : arange(dof_manager.getLocalSystemSize())) {
      if (not dof_manager.isLocalOrMasterDOF(d)) {
        idx.push_back(dof_manager.localToGlobalEquationNumber(d));
      }
    }
    AKANTU_DEBUG_ASSERT(
        idx.size() == (local_system_size - nb_local_dofs),
        "One developer that I will not name does not know how to count");

    PETSc_call(VecMPISetGhost, x, idx.size(), idx.data());
  }

  auto & is_ltog_mapping = dof_manager.getISLocalToGlobalMapping();

  PETSc_call(VecSetLocalToGlobalMapping, x, is_ltog_mapping);
}

/* -------------------------------------------------------------------------- */
void SolverVectorPETSc::clear() {
  PETSc_call(VecSet, x, 0.);
  applyModifications();
}

/* -------------------------------------------------------------------------- */
void SolverVectorPETSc::applyModifications() {
  PETSc_call(VecAssemblyBegin, x);
  PETSc_call(VecAssemblyEnd, x);
}

/* -------------------------------------------------------------------------- */
SolverVectorPETSc::operator const Array<Real> &()  const {
  const_cast<Array<Real> &>(cache).resize(local_size());
  
  auto xl = internal::make_petsc_local_vector(x);
  auto cachep = internal::make_petsc_wraped_vector(this->cache);

  PETSc_call(VecCopy, cachep, xl);
  return cache;
}

} // namespace akantu
