/**
 * @file   solver_vector_petsc.cc
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
#include "solver_vector_petsc.hh"
#include "dof_manager_petsc.hh"
#include "mpi_communicator_data.hh"
/* -------------------------------------------------------------------------- */
#include <numeric>
#include <petscvec.h>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
SolverVectorPETSc::SolverVectorPETSc(DOFManagerPETSc & dof_manager,
                                     const ID & id)
    : SolverVector(dof_manager, id), dof_manager(dof_manager) {
  auto && mpi_comm = dof_manager.getMPIComm();
  PETSc_call(VecCreate, mpi_comm, &x);
  detail::PETScSetName(x, id);

  PETSc_call(VecSetFromOptions, x);

  auto local_system_size = dof_manager.getLocalSystemSize();
  auto nb_local_dofs = dof_manager.getPureLocalSystemSize();
  PETSc_call(VecSetSizes, x, nb_local_dofs, PETSC_DECIDE);

  VecType vec_type;
  PETSc_call(VecGetType, x, &vec_type);
  if (std::string(vec_type) == std::string(VECMPI)) {
    PetscInt lowest_gidx, highest_gidx;
    PETSc_call(VecGetOwnershipRange, x, &lowest_gidx, &highest_gidx);

    std::vector<PetscInt> ghost_idx;
    for (auto && d : arange(local_system_size)) {
      int gidx = dof_manager.localToGlobalEquationNumber(d);
      if (gidx != -1) {
        if ((gidx < lowest_gidx) or (gidx >= highest_gidx)) {
          ghost_idx.push_back(gidx);
        }
      }
    }

    PETSc_call(VecMPISetGhost, x, ghost_idx.size(), ghost_idx.data());
  } else {
    std::vector<int> idx(nb_local_dofs);
    std::iota(idx.begin(), idx.end(), 0);
    ISLocalToGlobalMapping is;
    PETSc_call(ISLocalToGlobalMappingCreate, PETSC_COMM_SELF, 1, idx.size(),
               idx.data(), PETSC_COPY_VALUES, &is);
    PETSc_call(VecSetLocalToGlobalMapping, x, is);
    PETSc_call(ISLocalToGlobalMappingDestroy, &is);
  }
}

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
void SolverVectorPETSc::printself(std::ostream & stream, int indent) const {
  std::string space(indent, AKANTU_INDENT);
  stream << space << "SolverVectorPETSc [" << std::endl;
  stream << space << " + id: " << id << std::endl;
  PETSc_call(PetscViewerPushFormat, PETSC_VIEWER_STDOUT_WORLD,
             PETSC_VIEWER_ASCII_INDEX);
  PETSc_call(VecView, x, PETSC_VIEWER_STDOUT_WORLD);
  PETSc_call(PetscViewerPopFormat, PETSC_VIEWER_STDOUT_WORLD);
  stream << space << "]" << std::endl;
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
  // the arrays are destroyed and recreated in the dof manager
  // resize is so not implemented
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
  updateGhost();
}

/* -------------------------------------------------------------------------- */
void SolverVectorPETSc::updateGhost() {
  Vec x_ghosted{nullptr};
  PETSc_call(VecGhostGetLocalForm, x, &x_ghosted);
  if (x_ghosted) {
    PETSc_call(VecGhostUpdateBegin, x, INSERT_VALUES, SCATTER_FORWARD);
    PETSc_call(VecGhostUpdateEnd, x, INSERT_VALUES, SCATTER_FORWARD);
  }
  PETSc_call(VecGhostRestoreLocalForm, x, &x_ghosted);
}

/* -------------------------------------------------------------------------- */
void SolverVectorPETSc::getValues(const Array<Int> & idx,
                                  Array<Real> & values) const {
  if (idx.size() == 0)
    return;

  ISLocalToGlobalMapping is_ltog_map;
  PETSc_call(VecGetLocalToGlobalMapping, x, &is_ltog_map);

  PetscInt n;
  Array<PetscInt> lidx(idx.size());
  PETSc_call(ISGlobalToLocalMappingApply, is_ltog_map, IS_GTOLM_MASK,
             idx.size(), idx.storage(), &n, lidx.storage());

  getValuesLocal(lidx, values);
}
/* -------------------------------------------------------------------------- */
void SolverVectorPETSc::getValuesLocal(const Array<Int> & idx,
                                       Array<Real> & values) const {
  if (idx.size() == 0)
    return;

  Vec x_ghosted{nullptr};
  PETSc_call(VecGhostGetLocalForm, x, &x_ghosted);
  // VecScatterBegin(scatter, x, x_local, INSERT_VALUES, SCATTER_FORWARD);
  // VecScatterEnd(scatter, x, x_local, INSERT_VALUES, SCATTER_FORWARD);

  if (not x_ghosted) {
    const PetscScalar * array;
    PETSc_call(VecGetArrayRead, x, &array);

    for (auto && data : zip(idx, make_view(values))) {
      auto i = std::get<0>(data);
      if (i != -1) {
        std::get<1>(data) = array[i];
      }
    }

    PETSc_call(VecRestoreArrayRead, x, &array);
    return;
  }

  PETSc_call(VecSetOption, x_ghosted, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
  PETSc_call(VecGetValues, x_ghosted, idx.size(), idx.storage(),
             values.storage());
  PETSc_call(VecGhostRestoreLocalForm, x, &x_ghosted);
}

/* -------------------------------------------------------------------------- */
void SolverVectorPETSc::addValues(const Array<Int> & gidx,
                                  const Array<Real> & values,
                                  Real scale_factor) {
  Real * to_add = values.storage();
  Array<Real> scaled_array;
  if (scale_factor != 1.) {
    scaled_array.copy(values, false);
    scaled_array *= scale_factor;
    to_add = scaled_array.storage();
  }

  PETSc_call(VecSetOption, x, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
  PETSc_call(VecSetValues, x, gidx.size(), gidx.storage(), to_add, ADD_VALUES);

  applyModifications();
}

/* -------------------------------------------------------------------------- */
void SolverVectorPETSc::addValuesLocal(const Array<Int> & lidx,
                                       const Array<Real> & values,
                                       Real scale_factor) {
  Vec x_ghosted{nullptr};
  PETSc_call(VecGhostGetLocalForm, x, &x_ghosted);

  if (not x_ghosted) {
    Real * to_add = values.storage();
    Array<Real> scaled_array;
    if (scale_factor != 1.) {
      scaled_array.copy(values, false);
      scaled_array *= scale_factor;
      to_add = scaled_array.storage();
    }

    PETSc_call(VecSetOption, x, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
    PETSc_call(VecSetValuesLocal, x, lidx.size(), lidx.storage(), to_add,
               ADD_VALUES);
    return;
  }

  PETSc_call(VecGhostRestoreLocalForm, x, &x_ghosted);

  ISLocalToGlobalMapping is_ltog_map;
  PETSc_call(VecGetLocalToGlobalMapping, x, &is_ltog_map);

  Array<Int> gidx(lidx.size());
  PETSc_call(ISLocalToGlobalMappingApply, is_ltog_map, lidx.size(),
             lidx.storage(), gidx.storage());
  addValues(gidx, values, scale_factor);
}

/* -------------------------------------------------------------------------- */
SolverVectorPETSc::operator const Array<Real> &() const {
  const_cast<Array<Real> &>(cache).resize(local_size());

  auto xl = internal::make_petsc_local_vector(x);
  auto cachep = internal::make_petsc_wraped_vector(this->cache);

  PETSc_call(VecCopy, cachep, xl);
  return cache;
}

/* -------------------------------------------------------------------------- */
SolverVectorPETSc & SolverVectorPETSc::operator=(const SolverVectorPETSc & y) {
  if (size() != y.size()) {
    PETSc_call(VecDuplicate, y, &x);
  }

  PETSc_call(VecCopy, y.x, x);
  release_ = y.release_;
  return *this;
}

/* -------------------------------------------------------------------------- */
SolverVector & SolverVectorPETSc::operator=(const SolverVector & y) {
  const auto & y_ = aka::as_type<SolverVectorPETSc>(y);
  return operator=(y_);
}

/* -------------------------------------------------------------------------- */
SolverVector & SolverVectorPETSc::operator+(const SolverVector & y) {
  auto & y_ = aka::as_type<SolverVectorPETSc>(y);
  PETSc_call(VecAXPY, x, 1., y_.x);
  release_ = y_.release_;
  return *this;
}

} // namespace akantu
