/**
 * @file   dof_manager_default.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Aug 18 2015
 * @date last modification: Wed Jan 31 2018
 *
 * @brief  Default implementation of the dof manager
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
#include <functional>
#include <unordered_map>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_DOF_MANAGER_DEFAULT_HH__
#define __AKANTU_DOF_MANAGER_DEFAULT_HH__

namespace akantu {
class SparseMatrixAIJ;
class NonLinearSolverDefault;
class TimeStepSolverDefault;
class DOFSynchronizer;
class SolverVectorDefault;
} // namespace akantu

namespace akantu {

class DOFManagerDefault : public DOFManager {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  DOFManagerDefault(const ID & id = "dof_manager_default",
                    const MemoryID & memory_id = 0);
  DOFManagerDefault(Mesh & mesh, const ID & id = "dof_manager_default",
                    const MemoryID & memory_id = 0);
  ~DOFManagerDefault() override;

protected:
  struct DOFDataDefault : public DOFData {
    explicit DOFDataDefault(const ID & dof_id);
  };

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  // /// register an array of degree of freedom
  // void registerDOFs(const ID & dof_id, Array<Real> & dofs_array,
  //                   const DOFSupportType & support_type) override;

  // /// the dof as an implied type of _dst_nodal and is defined only on a
  // subset
  // /// of nodes
  // void registerDOFs(const ID & dof_id, Array<Real> & dofs_array,
  //                   const ID & group_support) override;

  /**
   * Assemble elementary values to the global matrix. The dof number is
   * implicitly considered as conn(el, n) * nb_nodes_per_element + d.
   * With 0 < n < nb_nodes_per_element and 0 < d < nb_dof_per_node
   **/
  void assembleElementalMatricesToMatrix(
      const ID & matrix_id, const ID & dof_id,
      const Array<Real> & elementary_mat, const ElementType & type,
      const GhostType & ghost_type, const MatrixType & elemental_matrix_type,
      const Array<UInt> & filter_elements) override;

  void assembleMatMulVectToGlobalArray(const ID & dof_id, const ID & A_id,
                                       const Array<Real> & x,
                                       Array<Real> & array,
                                       Real scale_factor = 1.);

  void assembleMatMulVectToArray(const ID & dof_id, const ID & A_id,
                                 const Array<Real> & x, Array<Real> & array,
                                 Real scale_factor = 1.) override;

  /// multiply a vector by a matrix and assemble the result to the residual
  void assembleMatMulVectToResidual(const ID & dof_id, const ID & A_id,
                                    const Array<Real> & x,
                                    Real scale_factor = 1) override;

  /// multiply a vector by a lumped matrix and assemble the result to the
  /// residual
  void assembleLumpedMatMulVectToResidual(const ID & dof_id, const ID & A_id,
                                          const Array<Real> & x,
                                          Real scale_factor = 1) override;

  /// assemble coupling terms between to dofs
  void assemblePreassembledMatrix(const ID & dof_id_m, const ID & dof_id_n,
                                  const ID & matrix_id,
                                  const TermsToAssemble & terms) override;

protected:
  void assembleToGlobalArray(const ID & dof_id,
                             const Array<Real> & array_to_assemble,
                             SolverVector & global_array,
                             Real scale_factor) override;

  template <typename T>
  void assembleToGlobalArray(const ID & dof_id,
                             const Array<T> & array_to_assemble,
                             Array<T> & global_array, T scale_factor);

  void getArrayPerDOFs(const ID & dof_id, const SolverVector & global,
                       Array<Real> & local) override;

  template <typename T>
  void getArrayPerDOFs(const ID & dof_id,
                                          const Array<T> & global_array,
                                          Array<T> & local_array) const;
  void makeConsistentForPeriodicity(const ID & dof_id,
                                    SolverVector & array) override;

public:
  /// update the global dofs vector
  virtual void updateGlobalBlockedDofs();

  /// apply boundary conditions to jacobian matrix
  void applyBoundary(const ID & matrix_id = "J") override;

private:
  /// Add a symmetric matrices to a symmetric sparse matrix
  void addSymmetricElementalMatrixToSymmetric(
      SparseMatrixAIJ & matrix, const Matrix<Real> & element_mat,
      const Vector<Int> & equation_numbers, UInt max_size);

  /// Add a unsymmetric matrices to a symmetric sparse matrix (i.e. cohesive
  /// elements)
  void addUnsymmetricElementalMatrixToSymmetric(
      SparseMatrixAIJ & matrix, const Matrix<Real> & element_mat,
      const Vector<Int> & equation_numbers, UInt max_size);

  /// Add a matrices to a unsymmetric sparse matrix
  void addElementalMatrixToUnsymmetric(SparseMatrixAIJ & matrix,
                                       const Matrix<Real> & element_mat,
                                       const Vector<Int> & equation_numbers,
                                       UInt max_size);

  void addToProfile(const ID & matrix_id, const ID & dof_id,
                    const ElementType & type, const GhostType & ghost_type);

  /* ------------------------------------------------------------------------ */
  /* MeshEventHandler interface                                               */
  /* ------------------------------------------------------------------------ */
protected:
  std::tuple<UInt, UInt, UInt>
  registerDOFsInternal(const ID & dof_id, Array<Real> & dofs_array) override;

  // std::pair<UInt, UInt>
  // updateNodalDOFs(const ID & dof_id, const Array<UInt> & nodes_list)
  // override;

  void resizeGlobalArrays() override;

public:
  /// function to implement to react on  akantu::NewNodesEvent
  void onNodesAdded(const Array<UInt> & nodes_list,
                    const NewNodesEvent & event) override;

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
  SparseMatrixAIJ & getMatrix(const ID & matrix_id);

  /// Get an instance of a new lumped matrix
  SolverVector & getNewLumpedMatrix(const ID & matrix_id) override;

  /* ------------------------------------------------------------------------ */
  /* Non Linear Solver                                                        */
  /* ------------------------------------------------------------------------ */
  /// Get instance of a non linear solver
  NonLinearSolver & getNewNonLinearSolver(
      const ID & nls_solver_id,
      const NonLinearSolverType & _non_linear_solver_type) override;

  /* ------------------------------------------------------------------------ */
  /* Time-Step Solver                                                         */
  /* ------------------------------------------------------------------------ */
  /// Get instance of a time step solver
  TimeStepSolver &
  getNewTimeStepSolver(const ID & id, const TimeStepSolverType & type,
                       NonLinearSolver & non_linear_solver) override;

  /* ------------------------------------------------------------------------ */
private:
  /// Get the solution array
  Array<Real> & getSolutionArray();

  /// Get the residual array
  const Array<Real> & getResidualArray() const;

  /// Get the residual array
  Array<Real> & getResidualArray();

public:
  /// Get the blocked dofs array
  AKANTU_GET_MACRO(GlobalBlockedDOFs, global_blocked_dofs, const Array<bool> &);
  /// Get the blocked dofs array
  AKANTU_GET_MACRO(PreviousGlobalBlockedDOFs, previous_global_blocked_dofs,
                   const Array<bool> &);

  /// Get the equation numbers corresponding to a dof_id. This might be used to
  /// access the matrix.
  inline const Array<Int> & getLocalEquationsNumbers(const ID & dof_id) const;

  /// get the array of dof types (use only if you know what you do...)
  inline const Array<UInt> & getDOFsAssociatedNodes(const ID & dof_id) const;

  /// access the internal dof_synchronizer
  AKANTU_GET_MACRO_NOT_CONST(Synchronizer, *synchronizer, DOFSynchronizer &);

  /// access the internal dof_synchronizer
  bool hasSynchronizer() const { return synchronizer != nullptr; }

protected:
  std::unique_ptr<DOFData> getNewDOFData(const ID & dof_id) override;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  using DOFToMatrixProfile =
      std::map<std::pair<ID, ID>,
               std::vector<std::pair<ElementType, GhostType>>>;

  /// contains the the dofs that where added to the profile of a given matrix.
  DOFToMatrixProfile matrix_profiled_dofs;

  /// rhs used only on root proc in case of parallel computing, this is the full
  /// gathered rhs array
  std::unique_ptr<Array<Real>> global_residual;

  /// blocked degree of freedom in the system equation corresponding to the
  /// different dofs
  Array<bool> global_blocked_dofs;

  /// blocked degree of freedom in the system equation corresponding to the
  /// different dofs
  Array<bool> previous_global_blocked_dofs;

  /// Memory cache, this is an array to keep the temporary memory needed for
  /// some operations, it is meant to be resized or cleared when needed
  Array<Real> data_cache;

  /// Release at last apply boundary on jacobian
  UInt jacobian_release{0};

  /// synchronizer to maintain coherency in dof fields
  std::unique_ptr<DOFSynchronizer> synchronizer;
};

} // namespace akantu

#include "dof_manager_default_inline_impl.cc"

#endif /* __AKANTU_DOF_MANAGER_DEFAULT_HH__ */
