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
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_DOF_MANAGER_PETSC_HH_
#define AKANTU_DOF_MANAGER_PETSC_HH_

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
namespace detail {
  template <typename T> void PETScSetName(T t, const ID & id) {
    PETSc_call(PetscObjectSetName, reinterpret_cast<PetscObject>(t),
               id.c_str());
  }
} // namespace detail
} // namespace akantu

namespace akantu {
class SparseMatrixPETSc;
class SolverVectorPETSc;
} // namespace akantu

namespace akantu {

class DOFManagerPETSc : public DOFManager {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  DOFManagerPETSc(const ID & id = "dof_manager_petsc");
  DOFManagerPETSc(Mesh & mesh, const ID & id = "dof_manager_petsc");

  ~DOFManagerPETSc() override = default;

protected:
  void init();

  struct DOFDataPETSc : public DOFData {
    explicit DOFDataPETSc(const ID & dof_id);

    /// petsc compressed version of local_equation_number
    Array<PetscInt> local_equation_number_petsc;

    Array<Int> & getLocalEquationsNumbers() override {
      return local_equation_number_petsc;
    }
  };

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void assembleToLumpedMatrix(const ID & /*dof_id*/,
                              Array<Real> & /*array_to_assemble*/,
                              const ID & /*lumped_mtx*/,
                              Real /*scale_factor*/ = 1.) override {
    AKANTU_TO_IMPLEMENT();
  }

  void assembleElementalMatricesToMatrix(
      const ID & /*matrix_id*/, const ID & /*dof_id*/,
      const Array<Real> & /*elementary_mat*/, ElementType /*type*/,
      GhostType /*ghost_type*/,
      const MatrixType & /*elemental_matrix_type*/,
      const Array<UInt> & /*filter_elements*/) override;

  void assembleMatMulVectToArray(const ID & /*dof_id*/, const ID & /*A_id*/,
                                 const Array<Real> & /*x*/,
                                 Array<Real> & /*array*/,
                                 Real /*scale_factor*/ = 1.) override;

  void assembleLumpedMatMulVectToResidual(const ID & /*dof_id*/,
                                          const ID & /*A_id*/,
                                          const Array<Real> & /*x*/,
                                          Real /*scale_factor*/ = 1) override {
    AKANTU_TO_IMPLEMENT();
  }

  void assemblePreassembledMatrix(const ID & /* dof_id_m*/,
                                  const ID & /*dof_id_n*/,
                                  const ID & /*matrix_id*/,
                                  const TermsToAssemble & /*terms*/) override;

protected:
  void assembleToGlobalArray(const ID & dof_id,
                             const Array<Real> & array_to_assemble,
                             SolverVector & global_array,
                             Real scale_factor) override;
  void getArrayPerDOFs(const ID & dof_id, const SolverVector & global,
                       Array<Real> & local) override;

  void makeConsistentForPeriodicity(const ID & dof_id,
                                    SolverVector & array) override;

  std::unique_ptr<DOFData> getNewDOFData(const ID & dof_id) override;

  std::tuple<UInt, UInt, UInt>
  registerDOFsInternal(const ID & dof_id, Array<Real> & dofs_array) override;

  void updateDOFsData(DOFDataPETSc & dof_data, UInt nb_new_local_dofs,
                      UInt nb_new_pure_local, UInt nb_node,
                      const std::function<UInt(UInt)> & getNode);

protected:
  void getLumpedMatrixPerDOFs(const ID & /*dof_id*/, const ID & /*lumped_mtx*/,
                              Array<Real> & /*lumped*/) override {}

  NonLinearSolver & getNewNonLinearSolver(
      const ID & nls_solver_id,
      const NonLinearSolverType & non_linear_solver_type) override;

  TimeStepSolver &
  getNewTimeStepSolver(const ID & id, const TimeStepSolverType & type,
                       NonLinearSolver & non_linear_solver,
                       SolverCallback & solver_callback) override;

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

  /// Get an instance of a new lumped matrix
  SolverVector & getNewLumpedMatrix(const ID & matrix_id) override;

  /// Get the blocked dofs array
  //  AKANTU_GET_MACRO(BlockedDOFs, blocked_dofs, const Array<bool> &);
  AKANTU_GET_MACRO(MPIComm, mpi_communicator, MPI_Comm);

  AKANTU_GET_MACRO_NOT_CONST(ISLocalToGlobalMapping, is_ltog_map,
                             ISLocalToGlobalMapping &);

  SolverVectorPETSc & getSolution();
  const SolverVectorPETSc & getSolution() const;

  SolverVectorPETSc & getResidual();
  const SolverVectorPETSc & getResidual() const;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  using PETScMatrixMap = std::map<ID, SparseMatrixPETSc *>;
  using PETScLumpedMatrixMap = std::map<ID, SolverVectorPETSc *>;

  /// list of matrices registered to the dof manager
  PETScMatrixMap petsc_matrices;

  /// list of lumped matrices registered
  PETScLumpedMatrixMap petsc_lumped_matrices;

  /// PETSc local to global mapping of dofs
  ISLocalToGlobalMapping is_ltog_map{nullptr};

  /// Communicator associated to PETSc
  MPI_Comm mpi_communicator;

  /// list of the dof ids to be able to always iterate in the same order
  std::vector<ID> dofs_ids;
};

/* -------------------------------------------------------------------------- */

} // namespace akantu

#endif /* AKANTU_DOF_MANAGER_PETSC_HH_ */
