/**
 * @file   dof_manager_petsc.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Aug 18 2015
 * @date last modification: Wed Jan 31 2018
 *
 * @brief  PETSc implementation of the dof manager
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
#include "dof_manager.hh"
/* -------------------------------------------------------------------------- */
#include <petscis.h>
#include <petscvec.h>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_DOF_MANAGER_PETSC_HH__
#define __AKANTU_DOF_MANAGER_PETSC_HH__

#define PETSc_call(func, ...)                                                  \
  do {                                                                         \
    auto ierr = func(__VA_ARGS__);                                             \
    if (PetscUnlikely(ierr != 0)) {                                            \
      const char * desc;                                                       \
      PetscErrorMessage(ierr, &desc, nullptr);                                 \
      AKANTU_EXCEPTION("Error in PETSc call to \'" << #func                    \
                                                   << "\': " << desc);         \
    }                                                                          \
  } while (false)

namespace akantu {
class SparseMatrixPETSc;
}

namespace akantu {

class DOFManagerPETSc : public DOFManager {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  DOFManagerPETSc(const ID & id = "dof_manager_petsc",
                  const MemoryID & memory_id = 0);
  DOFManagerPETSc(Mesh & mesh, const ID & id = "dof_manager_petsc",
                  const MemoryID & memory_id = 0);

  virtual ~DOFManagerPETSc();

protected:
  void init();

  struct DOFDataPETSc : public DOFData {
    explicit DOFDataPETSc(const ID & dof_id);

    /// local equation numbers in PETSc type
    IS is;
  };

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void registerDOFs(const ID & dof_id, Array<Real> & dofs_array,
                    DOFSupportType & support_type);

  void assembleToResidual(const ID & dof_id, Array<Real> & array_to_assemble,
                          Real scale_factor = 1.) override;

  void assembleToLumpedMatrix(const ID & /*dof_id*/,
                              Array<Real> & /*array_to_assemble*/,
                              const ID & /*lumped_mtx*/,
                              Real /*scale_factor*/ = 1.) {
    AKANTU_TO_IMPLEMENT();
  }

  void assembleElementalMatricesToMatrix(
      const ID & /*matrix_id*/, const ID & /*dof_id*/,
      const Array<Real> & /*elementary_mat*/, const ElementType & /*type*/,
      const GhostType & /*ghost_type*/,
      const MatrixType & /*elemental_matrix_type*/,
      const Array<UInt> & /*filter_elements*/) override {
    AKANTU_TO_IMPLEMENT();
  }

  void assembleMatMulVectToArray(const ID & /*dof_id*/, const ID & /*A_id*/,
                                 const Array<Real> & /*x*/,
                                 Array<Real> & /*array*/,
                                 Real /*scale_factor*/ = 1.) override {
    AKANTU_TO_IMPLEMENT();
  }

  void assembleMatMulVectToResidual(const ID & /*dof_id*/, const ID & /*A_id*/,
                                    const Array<Real> & /*x*/,
                                    Real /*scale_factor*/ = 1) override {
    AKANTU_TO_IMPLEMENT();
  }

  void assembleLumpedMatMulVectToResidual(const ID & /*dof_id*/,
                                          const ID & /*A_id*/,
                                          const Array<Real> & /*x*/,
                                          Real /*scale_factor*/ = 1) override {
    AKANTU_TO_IMPLEMENT();
  }

  void assemblePreassembledMatrix(const ID & /* dof_id_m*/,
                                  const ID & /*dof_id_n*/,
                                  const ID & /*matrix_id*/,
                                  const TermsToAssemble & /*terms*/) override {
    AKANTU_TO_IMPLEMENT();
  }

  void clearResidual() override;
  void clearMatrix(const ID & mtx) override;
  void clearLumpedMatrix(const ID & mtx) override;

  void applyBoundary(const ID & /*matrix_id*/ = "J") override {
    AKANTU_TO_IMPLEMENT();
  }

protected:
  void getLumpedMatrixPerDOFs(const ID & dof_id, const ID & lumped_mtx,
                              Array<Real> & lumped) override;

  void getSolutionPerDOFs(const ID & dof_id,
                          Array<Real> & solution_array) override;

  NonLinearSolver & getNewNonLinearSolver(
      const ID & nls_solver_id,
      const NonLinearSolverType & non_linear_solver_type) override;

  TimeStepSolver &
  getNewTimeStepSolver(const ID & time_step_solver_id,
                       const TimeStepSolverType & type,
                       NonLinearSolver & non_linear_solver) override;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// Get an instance of a new SparseMatrix
  SparseMatrix & getNewMatrix(const ID & matrix_id,
                              const MatrixType & matrix_type) override;

  /// Get an instance of a new SparseMatrix as a copy of the SparseMatrix
  /// matrix_to_copy_id
  SparseMatrix & getNewMatrix(const ID & matrix_id,
                              const ID & matrix_to_copy_id) override;

  /// Get the reference of an existing matrix
  SparseMatrixPETSc & getMatrix(const ID & matrix_id);

  /// Get the reference of an existing matrix
  Vec & getLumpedMatrix(const ID & /*matrix_id*/) {
    AKANTU_TO_IMPLEMENT();
  }

  /// Get the solution array
  AKANTU_GET_MACRO_NOT_CONST(GlobalSolution, this->solution, Vec &);
  /// Get the residual array
  AKANTU_GET_MACRO_NOT_CONST(Residual, this->residual, Vec &);

  /// Get the blocked dofs array
  //  AKANTU_GET_MACRO(BlockedDOFs, blocked_dofs, const Array<bool> &);
  AKANTU_GET_MACRO(MPIComm, mpi_communicator, MPI_Comm);

  AKANTU_GET_MACRO_NOT_CONST(ISLocalToGlobalMapping, is_ltog_map,
                             ISLocalToGlobalMapping &);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  using PETScMatrixMap = std::map<ID, SparseMatrixPETSc *>;
  using PETScLumpedMatrixMap = std::map<ID, Vec>;

  /// list of matrices registered to the dof manager
  PETScMatrixMap petsc_matrices;

  /// list of lumped matrices registered
  PETScLumpedMatrixMap petsc_lumped_matrices;

  /// PETSc version of the solution
  Vec solution;

  /// PETSc version of the residual
  Vec residual;

  /// PETSc local to global mapping of dofs
  ISLocalToGlobalMapping is_ltog_map;

  /// Communicator associated to PETSc
  MPI_Comm mpi_communicator;

  /// counter of instances to know when to finalize
  static int petsc_dof_manager_instances;
};

} // namespace akantu

#endif /* __AKANTU_DOF_MANAGER_PETSC_HH__ */
