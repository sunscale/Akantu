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

    /// associated node for _dst_nodal dofs only
    Array<UInt> associated_nodes;
  };

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
private:
  void registerDOFsInternal(const ID & dof_id, UInt nb_dofs,
                            UInt nb_pure_local_dofs);

public:
  /// register an array of degree of freedom
  void registerDOFs(const ID & dof_id, Array<Real> & dofs_array,
                    const DOFSupportType & support_type) override;

  /// the dof as an implied type of _dst_nodal and is defined only on a subset
  /// of nodes
  void registerDOFs(const ID & dof_id, Array<Real> & dofs_array,
                    const ID & group_support) override;

  /// Assemble an array to the global residual array
  void assembleToResidual(const ID & dof_id,
                          const Array<Real> & array_to_assemble,
                          Real scale_factor = 1.) override;

  /// Assemble an array to the global lumped matrix array
  void assembleToLumpedMatrix(const ID & dof_id,
                              const Array<Real> & array_to_assemble,
                              const ID & lumped_mtx,
                              Real scale_factor = 1.) override;
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
  /// Assemble an array to the global residual array
  template <typename T>
  void assembleToGlobalArray(const ID & dof_id,
                             const Array<T> & array_to_assemble,
                             Array<T> & global_array, T scale_factor);

public:
  /// clear the residual
  void clearResidual() override;

  /// sets the matrix to 0
  void clearMatrix(const ID & mtx) override;

  /// sets the lumped matrix to 0
  void clearLumpedMatrix(const ID & mtx) override;

  /// update the global dofs vector
  virtual void updateGlobalBlockedDofs();

  /// apply boundary conditions to jacobian matrix
  virtual void applyBoundary(const ID & matrix_id = "J");

  // void getEquationsNumbers(const ID & dof_id,
  //                          Array<UInt> & equation_numbers) override;

protected:
  /// Get local part of an array corresponding to a given dofdata
  template <typename T>
  void getArrayPerDOFs(const ID & dof_id, const Array<T> & global_array,
                       Array<T> & local_array) const;

  /// Get the part of the solution corresponding to the dof_id
  void getSolutionPerDOFs(const ID & dof_id,
                          Array<Real> & solution_array) override;

  /// fill a Vector with the equation numbers corresponding to the given
  /// connectivity
  inline void extractElementEquationNumber(
      const Array<UInt> & equation_numbers, const Vector<UInt> & connectivity,
      UInt nb_degree_of_freedom, Vector<UInt> & local_equation_number);

public:
  /// extract a lumped matrix part corresponding to a given dof
  void getLumpedMatrixPerDOFs(const ID & dof_id, const ID & lumped_mtx,
                              Array<Real> & lumped) override;

private:
  /// Add a symmetric matrices to a symmetric sparse matrix
  void addSymmetricElementalMatrixToSymmetric(
      SparseMatrixAIJ & matrix, const Matrix<Real> & element_mat,
      const Vector<UInt> & equation_numbers, UInt max_size);

  /// Add a unsymmetric matrices to a symmetric sparse matrix (i.e. cohesive
  /// elements)
  void addUnsymmetricElementalMatrixToSymmetric(
      SparseMatrixAIJ & matrix, const Matrix<Real> & element_mat,
      const Vector<UInt> & equation_numbers, UInt max_size);

  /// Add a matrices to a unsymmetric sparse matrix
  void addElementalMatrixToUnsymmetric(SparseMatrixAIJ & matrix,
                                       const Matrix<Real> & element_mat,
                                       const Vector<UInt> & equation_numbers,
                                       UInt max_size);

  void addToProfile(const ID & matrix_id, const ID & dof_id,
                    const ElementType & type, const GhostType & ghost_type);

  /* ------------------------------------------------------------------------ */
  /* MeshEventHandler interface                                               */
  /* ------------------------------------------------------------------------ */
protected:
  std::pair<UInt, UInt>
  updateNodalDOFs(const ID & dof_id, const Array<UInt> & nodes_list) override;

private:
  void updateDOFsData(DOFDataDefault & dof_data, UInt nb_new_local_dofs,
                      UInt nb_new_pure_local,
                      const std::function<UInt(UInt)> & getNode);

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
  /// Get the solution array
  AKANTU_GET_MACRO_NOT_CONST(GlobalSolution, global_solution, Array<Real> &);
  /// Set the global solution array
  void setGlobalSolution(const Array<Real> & solution);

  /// Get the global residual array across processors (SMP call)
  const Array<Real> & getGlobalResidual();

  /// Get the residual array
  const Array<Real> & getResidual() const;

  /// Get the blocked dofs array
  AKANTU_GET_MACRO(GlobalBlockedDOFs, global_blocked_dofs, const Array<bool> &);
  /// Get the blocked dofs array
  AKANTU_GET_MACRO(PreviousGlobalBlockedDOFs, previous_global_blocked_dofs,
                   const Array<bool> &);
  /// Get the location type of a given dof
  inline bool isLocalOrMasterDOF(UInt local_dof_num);

  /// Answer to the question is a dof a slave dof ?
  inline bool isSlaveDOF(UInt local_dof_num);

  /// get the equation numbers (in local numbering) corresponding to a dof ID
  inline const Array<UInt> & getLocalEquationNumbers(const ID & dof_id) const;

  /// tells if the dof manager knows about a global dof
  bool hasGlobalEquationNumber(UInt global) const;

  /// return the local index of the global equation number
  inline UInt globalToLocalEquationNumber(UInt global) const;

  /// converts local equation numbers to global equation numbers;
  inline UInt localToGlobalEquationNumber(UInt local) const;

  /// get the array of dof types (use only if you know what you do...)
  inline Int getDOFType(UInt local_id) const;

  /// get the array of dof types (use only if you know what you do...)
  inline const Array<UInt> & getDOFsAssociatedNodes(const ID & dof_id) const;

  /// access the internal dof_synchronizer
  AKANTU_GET_MACRO_NOT_CONST(Synchronizer, *synchronizer, DOFSynchronizer &);

  /// access the internal dof_synchronizer
  bool hasSynchronizer() const { return synchronizer != nullptr; }

protected:
  DOFData & getNewDOFData(const ID & dof_id) override;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  // using AIJMatrixMap = std::map<ID, std::unique_ptr<SparseMatrixAIJ>>;
  // using DefaultNonLinearSolversMap =
  //     std::map<ID, std::unique_ptr<NonLinearSolverDefault>>;
  // using DefaultTimeStepSolversMap =
  //     std::map<ID, std::unique_ptr<TimeStepSolverDefault>>;

  using DOFToMatrixProfile =
      std::map<std::pair<ID, ID>,
               std::vector<std::pair<ElementType, GhostType>>>;

  /// contains the the dofs that where added to the profile of a given matrix.
  DOFToMatrixProfile matrix_profiled_dofs;

  /// rhs to the system of equation corresponding to the residual linked to the
  /// different dofs
  Array<Real> residual;

  /// rhs used only on root proc in case of parallel computing, this is the full
  /// gathered rhs array
  std::unique_ptr<Array<Real>> global_residual;

  /// solution of the system of equation corresponding to the different dofs
  Array<Real> global_solution;

  /// blocked degree of freedom in the system equation corresponding to the
  /// different dofs
  Array<bool> global_blocked_dofs;

  /// blocked degree of freedom in the system equation corresponding to the
  /// different dofs
  Array<bool> previous_global_blocked_dofs;

  /// define the dofs type, local, shared, ghost
  Array<Int> dofs_type;

  /// Map of the different matrices stored in the dof manager
  // AIJMatrixMap aij_matrices;

  /// Map of the different time step solvers stored with there real type
  // DefaultTimeStepSolversMap default_time_step_solver_map;

  /// Memory cache, this is an array to keep the temporary memory needed for
  /// some operations, it is meant to be resized or cleared when needed
  Array<Real> data_cache;

  /// Release at last apply boundary on jacobian
  UInt jacobian_release{0};

  /// equation number in global numbering
  Array<UInt> global_equation_number;

  using equation_numbers_map = std::unordered_map<UInt, UInt>;

  /// dual information of global_equation_number
  equation_numbers_map global_to_local_mapping;

  /// accumulator to know what would be the next global id to use
  UInt first_global_dof_id{0};

  /// synchronizer to maintain coherency in dof fields
  std::unique_ptr<DOFSynchronizer> synchronizer;
};

} // namespace akantu

#include "dof_manager_default_inline_impl.cc"

#endif /* __AKANTU_DOF_MANAGER_DEFAULT_HH__ */
