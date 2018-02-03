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

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// register an array of degree of freedom
  void registerDOFs(const ID & dof_id, Array<Real> & dofs_array,
                    DOFSupportType & support_type);

  /// Assemble an array to the global residual array
  virtual void assembleToResidual(const ID & dof_id,
                                  const Array<Real> & array_to_assemble,
                                  Real scale_factor = 1.);

  /**
   * Assemble elementary values to the global residual array. The dof number is
   * implicitly considered as conn(el, n) * nb_nodes_per_element + d.
   * With 0 < n < nb_nodes_per_element and 0 < d < nb_dof_per_node
   **/
  virtual void assembleElementalArrayResidual(
      const ID & dof_id, const Array<Real> & array_to_assemble,
      const ElementType & type, const GhostType & ghost_type,
      Real scale_factor = 1.);
  /**
   * Assemble elementary values to the global residual array. The dof number is
   * implicitly considered as conn(el, n) * nb_nodes_per_element + d.
   * With 0 < n < nb_nodes_per_element and 0 < d < nb_dof_per_node
   **/
  virtual void assembleElementalMatricesToMatrix(
      const ID & matrix_id, const ID & dof_id,
      const Array<Real> & elementary_mat, const ElementType & type,
      const GhostType & ghost_type, const MatrixType & elemental_matrix_type,
      const Array<UInt> & filter_elements);

protected:
  /// Get the part of the solution corresponding to the dof_id
  virtual void getSolutionPerDOFs(const ID & dof_id,
                                  Array<Real> & solution_array);

private:
  /// Add a symmetric matrices to a symmetric sparse matrix
  inline void addSymmetricElementalMatrixToSymmetric(
      SparseMatrixAIJ & matrix, const Matrix<Real> & element_mat,
      const Vector<UInt> & equation_numbers, UInt max_size);

  /// Add a unsymmetric matrices to a symmetric sparse matrix (i.e. cohesive
  /// elements)
  inline void addUnsymmetricElementalMatrixToSymmetric(
      SparseMatrixAIJ & matrix, const Matrix<Real> & element_mat,
      const Vector<UInt> & equation_numbers, UInt max_size);

  /// Add a matrices to a unsymmetric sparse matrix
  inline void addElementalMatrixToUnsymmetric(
      SparseMatrixAIJ & matrix, const Matrix<Real> & element_mat,
      const Vector<UInt> & equation_numbers, UInt max_size);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// Get an instance of a new SparseMatrix
  virtual SparseMatrix & getNewMatrix(const ID & matrix_id,
                                      const MatrixType & matrix_type);

  /// Get an instance of a new SparseMatrix as a copy of the SparseMatrix
  /// matrix_to_copy_id
  virtual SparseMatrix & getNewMatrix(const ID & matrix_id,
                                      const ID & matrix_to_copy_id);

  /// Get the reference of an existing matrix
  SparseMatrixPETSc & getMatrix(const ID & matrix_id);

  /// Get the solution array
  AKANTU_GET_MACRO_NOT_CONST(GlobalSolution, this->solution, Vec &);
  /// Get the residual array
  AKANTU_GET_MACRO_NOT_CONST(Residual, this->residual, Vec &);

  /// Get the blocked dofs array
  //  AKANTU_GET_MACRO(BlockedDOFs, blocked_dofs, const Array<bool> &);
  AKANTU_GET_MACRO(MPIComm, mpi_communicator, MPI_Comm);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  typedef std::map<ID, SparseMatrixPETSc *> PETScMatrixMap;

  /// list of matrices registered to the dof manager
  PETScMatrixMap petsc_matrices;

  /// PETSc version of the solution
  Vec solution;

  /// PETSc version of the residual
  Vec residual;

  /// PETSc local to global mapping of dofs
  ISLocalToGlobalMapping is_ltog;

  /// Communicator associated to PETSc
  MPI_Comm mpi_communicator;

  static int petsc_dof_manager_instances{0};
};

} // akantu

#endif /* __AKANTU_DOF_MANAGER_PETSC_HH__ */
