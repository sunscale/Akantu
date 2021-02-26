/**
 * @file   dof_manager.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Aug 18 2015
 * @date last modification: Wed Feb 21 2018
 *
 * @brief  Class handling the different types of dofs
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
#include "aka_factory.hh"
#include "aka_memory.hh"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */
#include <map>
#include <set>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_DOF_MANAGER_HH_
#define AKANTU_DOF_MANAGER_HH_

namespace akantu {
class TermsToAssemble;
class NonLinearSolver;
class TimeStepSolver;
class SparseMatrix;
class SolverVector;
class SolverCallback;
} // namespace akantu

namespace akantu {

class DOFManager : protected Memory, protected MeshEventHandler {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
protected:
  struct DOFData;

public:
  DOFManager(const ID & id = "dof_manager", const MemoryID & memory_id = 0);
  DOFManager(Mesh & mesh, const ID & id = "dof_manager",
             const MemoryID & memory_id = 0);
  ~DOFManager() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// register an array of degree of freedom
  virtual void registerDOFs(const ID & dof_id, Array<Real> & dofs_array,
                            const DOFSupportType & support_type);

  /// the dof as an implied type of _dst_nodal and is defined only on a subset
  /// of nodes
  virtual void registerDOFs(const ID & dof_id, Array<Real> & dofs_array,
                            const ID & support_group);

  /// register an array of previous values of the degree of freedom
  virtual void registerDOFsPrevious(const ID & dof_id,
                                    Array<Real> & dofs_array);

  /// register an array of increment of degree of freedom
  virtual void registerDOFsIncrement(const ID & dof_id,
                                     Array<Real> & dofs_array);

  /// register an array of derivatives for a particular dof array
  virtual void registerDOFsDerivative(const ID & dof_id, UInt order,
                                      Array<Real> & dofs_derivative);

  /// register array representing the blocked degree of freedoms
  virtual void registerBlockedDOFs(const ID & dof_id,
                                   Array<bool> & blocked_dofs);

  /// Assemble an array to the global residual array
  virtual void assembleToResidual(const ID & dof_id,
                                  Array<Real> & array_to_assemble,
                                  Real scale_factor = 1.);

  /// Assemble an array to the global lumped matrix array
  virtual void assembleToLumpedMatrix(const ID & dof_id,
                                      Array<Real> & array_to_assemble,
                                      const ID & lumped_mtx,
                                      Real scale_factor = 1.);

  /**
   * Assemble elementary values to a local array of the size nb_nodes *
   * nb_dof_per_node. The dof number is implicitly considered as
   * conn(el, n) * nb_nodes_per_element + d.
   * With 0 < n < nb_nodes_per_element and 0 < d < nb_dof_per_node
   **/
  virtual void assembleElementalArrayLocalArray(
      const Array<Real> & elementary_vect, Array<Real> & array_assembeled,
      ElementType type, GhostType ghost_type,
      Real scale_factor = 1.,
      const Array<UInt> & filter_elements = empty_filter);

  /**
   * Assemble elementary values to the global residual array. The dof number is
   * implicitly considered as conn(el, n) * nb_nodes_per_element + d.
   * With 0 < n < nb_nodes_per_element and 0 < d < nb_dof_per_node
   **/
  virtual void assembleElementalArrayToResidual(
      const ID & dof_id, const Array<Real> & elementary_vect,
      ElementType type, GhostType ghost_type,
      Real scale_factor = 1.,
      const Array<UInt> & filter_elements = empty_filter);

  /**
   * Assemble elementary values to a global array corresponding to a lumped
   * matrix
   */
  virtual void assembleElementalArrayToLumpedMatrix(
      const ID & dof_id, const Array<Real> & elementary_vect,
      const ID & lumped_mtx, ElementType type,
      GhostType ghost_type, Real scale_factor = 1.,
      const Array<UInt> & filter_elements = empty_filter);

  /**
   * Assemble elementary values to the global residual array. The dof number is
   * implicitly considered as conn(el, n) * nb_nodes_per_element + d.  With 0 <
   * n < nb_nodes_per_element and 0 < d < nb_dof_per_node
   **/
  virtual void assembleElementalMatricesToMatrix(
      const ID & matrix_id, const ID & dof_id,
      const Array<Real> & elementary_mat, ElementType type,
      GhostType ghost_type = _not_ghost,
      const MatrixType & elemental_matrix_type = _symmetric,
      const Array<UInt> & filter_elements = empty_filter) = 0;

  /// multiply a vector by a matrix and assemble the result to the residual
  virtual void assembleMatMulVectToArray(const ID & dof_id, const ID & A_id,
                                         const Array<Real> & x,
                                         Array<Real> & array,
                                         Real scale_factor = 1) = 0;

  /// multiply a vector by a lumped matrix and assemble the result to the
  /// residual
  virtual void assembleLumpedMatMulVectToResidual(const ID & dof_id,
                                                  const ID & A_id,
                                                  const Array<Real> & x,
                                                  Real scale_factor = 1) = 0;

  /// assemble coupling terms between to dofs
  virtual void assemblePreassembledMatrix(const ID & dof_id_m,
                                          const ID & dof_id_n,
                                          const ID & matrix_id,
                                          const TermsToAssemble & terms) = 0;

  /// multiply a vector by a matrix and assemble the result to the residual
  virtual void assembleMatMulVectToResidual(const ID & dof_id, const ID & A_id,
                                            const Array<Real> & x,
                                            Real scale_factor = 1);

  /// multiply the dofs by a matrix and assemble the result to the residual
  virtual void assembleMatMulDOFsToResidual(const ID & A_id,
                                            Real scale_factor = 1);

  /// updates the global blocked_dofs array
  virtual void updateGlobalBlockedDofs();

  /// sets the residual to 0
  virtual void zeroResidual();
  /// sets the matrix to 0
  virtual void zeroMatrix(const ID & mtx);
  /// sets the lumped matrix to 0
  virtual void zeroLumpedMatrix(const ID & mtx);

  virtual void applyBoundary(const ID & matrix_id = "J");
  // virtual void applyBoundaryLumped(const ID & matrix_id = "J");

  /// extract a lumped matrix part corresponding to a given dof
  virtual void getLumpedMatrixPerDOFs(const ID & dof_id, const ID & lumped_mtx,
                                      Array<Real> & lumped);

  /// splits the solution storage from a global view to the per dof storages
  void splitSolutionPerDOFs();

private:
  /// dispatch the creation of the dof data and register it
  DOFData & getNewDOFDataInternal(const ID & dof_id);

protected:
  /// common function to help registering dofs the return values are the add new
  /// numbers of local dofs, pure local dofs, and system size
  virtual std::tuple<UInt, UInt, UInt>
  registerDOFsInternal(const ID & dof_id, Array<Real> & dofs_array);

  /// minimum functionality to implement per derived version of the DOFManager
  /// to allow the splitSolutionPerDOFs function to work
  virtual void getSolutionPerDOFs(const ID & dof_id,
                                  Array<Real> & solution_array);

  /// fill a Vector with the equation numbers corresponding to the given
  /// connectivity
  static inline void extractElementEquationNumber(
      const Array<Int> & equation_numbers, const Vector<UInt> & connectivity,
      UInt nb_degree_of_freedom, Vector<Int> & element_equation_number);

  /// Assemble a array to a global one
  void assembleMatMulVectToGlobalArray(const ID & dof_id, const ID & A_id,
                                       const Array<Real> & x,
                                       SolverVector & array,
                                       Real scale_factor = 1.);

  /// common function that can be called by derived class with proper matrice
  /// types
  template <typename Mat>
  void assemblePreassembledMatrix_(Mat & A, const ID & dof_id_m,
                                   const ID & dof_id_n,
                                   const TermsToAssemble & terms);

  template <typename Mat>
  void assembleElementalMatricesToMatrix_(
      Mat & A, const ID & dof_id, const Array<Real> & elementary_mat,
      ElementType type, GhostType ghost_type,
      const MatrixType & elemental_matrix_type,
      const Array<UInt> & filter_elements);

  template <typename Vec>
  void assembleMatMulVectToArray_(const ID & dof_id, const ID & A_id,
                                  const Array<Real> & x, Array<Real> & array,
                                  Real scale_factor);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// Get the location type of a given dof
  inline bool isLocalOrMasterDOF(UInt local_dof_num);

  /// Answer to the question is a dof a slave dof ?
  inline bool isSlaveDOF(UInt local_dof_num);

  /// Answer to the question is a dof a slave dof ?
  inline bool isPureGhostDOF(UInt local_dof_num);

  /// tells if the dof manager knows about a global dof
  bool hasGlobalEquationNumber(Int global) const;

  /// return the local index of the global equation number
  inline Int globalToLocalEquationNumber(Int global) const;

  /// converts local equation numbers to global equation numbers;
  inline Int localToGlobalEquationNumber(Int local) const;

  /// get the array of dof types (use only if you know what you do...)
  inline NodeFlag getDOFFlag(Int local_id) const;

  /// Global number of dofs
  AKANTU_GET_MACRO(SystemSize, this->system_size, UInt);

  /// Local number of dofs
  AKANTU_GET_MACRO(LocalSystemSize, this->local_system_size, UInt);

  /// Pure local number of dofs
  AKANTU_GET_MACRO(PureLocalSystemSize, this->pure_local_system_size, UInt);

  /// Retrieve all the registered DOFs
  std::vector<ID> getDOFIDs() const;

  /* ------------------------------------------------------------------------ */
  /* DOFs and derivatives accessors                                          */
  /* ------------------------------------------------------------------------ */
  /// Get a reference to the registered dof array for a given id
  inline Array<Real> & getDOFs(const ID & dofs_id);

  /// Get the support type of a given dof
  inline DOFSupportType getSupportType(const ID & dofs_id) const;

  /// are the dofs registered
  inline bool hasDOFs(const ID & dof_id) const;

  /// Get a reference to the registered dof derivatives array for a given id
  inline Array<Real> & getDOFsDerivatives(const ID & dofs_id, UInt order);

  /// Does the dof has derivatives
  inline bool hasDOFsDerivatives(const ID & dofs_id, UInt order) const;

  /// Get a reference to the blocked dofs array registered for the given id
  inline const Array<bool> & getBlockedDOFs(const ID & dofs_id) const;

  /// Does the dof has a blocked array
  inline bool hasBlockedDOFs(const ID & dofs_id) const;

  /// Get a reference to the registered dof increment array for a given id
  inline Array<Real> & getDOFsIncrement(const ID & dofs_id);

  /// Does the dof has a increment array
  inline bool hasDOFsIncrement(const ID & dofs_id) const;

  /// Does the dof has a previous array
  inline Array<Real> & getPreviousDOFs(const ID & dofs_id);

  /// Get a reference to the registered dof array for previous step values a
  /// given id
  inline bool hasPreviousDOFs(const ID & dofs_id) const;

  /// saves the values from dofs to previous dofs
  virtual void savePreviousDOFs(const ID & dofs_id);

  /// Get a reference to the solution array registered for the given id
  inline const Array<Real> & getSolution(const ID & dofs_id) const;

  /// Get a reference to the solution array registered for the given id
  inline Array<Real> & getSolution(const ID & dofs_id);

  /// Get the blocked dofs array
  AKANTU_GET_MACRO(GlobalBlockedDOFs, global_blocked_dofs, const Array<Int> &);
  /// Get the blocked dofs array
  AKANTU_GET_MACRO(PreviousGlobalBlockedDOFs, previous_global_blocked_dofs,
                   const Array<Int> &);

  /* ------------------------------------------------------------------------ */
  /* Matrices accessors                                                       */
  /* ------------------------------------------------------------------------ */
  /// Get an instance of a new SparseMatrix
  virtual SparseMatrix & getNewMatrix(const ID & matrix_id,
                                      const MatrixType & matrix_type) = 0;

  /// Get an instance of a new SparseMatrix as a copy of the SparseMatrix
  /// matrix_to_copy_id
  virtual SparseMatrix & getNewMatrix(const ID & matrix_id,
                                      const ID & matrix_to_copy_id) = 0;

  /// Get the equation numbers corresponding to a dof_id. This might be used to
  /// access the matrix.
  inline const Array<Int> & getLocalEquationsNumbers(const ID & dof_id) const;

protected:
  /// get the array of dof types (use only if you know what you do...)
  inline const Array<UInt> & getDOFsAssociatedNodes(const ID & dof_id) const;

protected:
  /* ------------------------------------------------------------------------ */
  /// register a matrix
  SparseMatrix & registerSparseMatrix(const ID & matrix_id,
                                      std::unique_ptr<SparseMatrix> & matrix);

  /// register a lumped matrix (aka a Vector)
  SolverVector & registerLumpedMatrix(const ID & matrix_id,
                                      std::unique_ptr<SolverVector> & matrix);

  /// register a non linear solver instantiated by a derived class
  NonLinearSolver &
  registerNonLinearSolver(const ID & non_linear_solver_id,
                          std::unique_ptr<NonLinearSolver> & non_linear_solver);

  /// register a time step solver instantiated by a derived class
  TimeStepSolver &
  registerTimeStepSolver(const ID & time_step_solver_id,
                         std::unique_ptr<TimeStepSolver> & time_step_solver);

  template <class NLSType, class DMType>
  NonLinearSolver & registerNonLinearSolver(DMType & dm, const ID & id,
                                            const NonLinearSolverType & type) {
    ID non_linear_solver_id = this->id + ":nls:" + id;
    std::unique_ptr<NonLinearSolver> nls = std::make_unique<NLSType>(
        dm, type, non_linear_solver_id, this->memory_id);
    return this->registerNonLinearSolver(non_linear_solver_id, nls);
  }

  template <class TSSType, class DMType>
  TimeStepSolver & registerTimeStepSolver(DMType & dm, const ID & id,
                                          const TimeStepSolverType & type,
                                          NonLinearSolver & non_linear_solver,
                                          SolverCallback & solver_callback) {
    ID time_step_solver_id = this->id + ":tss:" + id;
    std::unique_ptr<TimeStepSolver> tss =
        std::make_unique<TSSType>(dm, type, non_linear_solver, solver_callback,
                                  time_step_solver_id, this->memory_id);
    return this->registerTimeStepSolver(time_step_solver_id, tss);
  }

  template <class MatType, class DMType>
  SparseMatrix & registerSparseMatrix(DMType & dm, const ID & id,
                                      const MatrixType & matrix_type) {
    ID matrix_id = this->id + ":mtx:" + id;
    std::unique_ptr<SparseMatrix> sm =
        std::make_unique<MatType>(dm, matrix_type, matrix_id);
    return this->registerSparseMatrix(matrix_id, sm);
  }

  template <class MatType>
  SparseMatrix & registerSparseMatrix(const ID & id,
                                      const ID & matrix_to_copy_id) {
    ID matrix_id = this->id + ":mtx:" + id;
    auto & sm_to_copy =
        aka::as_type<MatType>(this->getMatrix(matrix_to_copy_id));
    std::unique_ptr<SparseMatrix> sm =
        std::make_unique<MatType>(sm_to_copy, matrix_id);
    return this->registerSparseMatrix(matrix_id, sm);
  }

  template <class MatType, class DMType>
  SolverVector & registerLumpedMatrix(DMType & dm, const ID & id) {
    ID matrix_id = this->id + ":lumped_mtx:" + id;
    std::unique_ptr<SolverVector> sm = std::make_unique<MatType>(dm, matrix_id);
    return this->registerLumpedMatrix(matrix_id, sm);
  }

protected:
  virtual void makeConsistentForPeriodicity(const ID & dof_id,
                                            SolverVector & array) = 0;

  virtual void assembleToGlobalArray(const ID & dof_id,
                                     const Array<Real> & array_to_assemble,
                                     SolverVector & global_array,
                                     Real scale_factor) = 0;

public:
  /// extract degrees of freedom (identified by ID) from a global solver array
  virtual void getArrayPerDOFs(const ID & dof_id, const SolverVector & global,
                               Array<Real> & local) = 0;

  /// Get the reference of an existing matrix
  SparseMatrix & getMatrix(const ID & matrix_id);

  /// check if the given matrix exists
  bool hasMatrix(const ID & matrix_id) const;

  /// Get an instance of a new lumped matrix
  virtual SolverVector & getNewLumpedMatrix(const ID & matrix_id) = 0;
  /// Get the lumped version of a given matrix
  const SolverVector & getLumpedMatrix(const ID & matrix_id) const;
  /// Get the lumped version of a given matrix
  SolverVector & getLumpedMatrix(const ID & matrix_id);

  /// check if the given matrix exists
  bool hasLumpedMatrix(const ID & matrix_id) const;

  /* ------------------------------------------------------------------------ */
  /* Non linear system solver                                                 */
  /* ------------------------------------------------------------------------ */
  /// Get instance of a non linear solver
  virtual NonLinearSolver & getNewNonLinearSolver(
      const ID & nls_solver_id,
      const NonLinearSolverType & _non_linear_solver_type) = 0;

  /// get instance of a non linear solver
  virtual NonLinearSolver & getNonLinearSolver(const ID & nls_solver_id);

  /// check if the given solver exists
  bool hasNonLinearSolver(const ID & solver_id) const;

  /* ------------------------------------------------------------------------ */
  /* Time-Step Solver                                                         */
  /* ------------------------------------------------------------------------ */
  /// Get instance of a time step solver
  virtual TimeStepSolver &
  getNewTimeStepSolver(const ID & time_step_solver_id,
                       const TimeStepSolverType & type,
                       NonLinearSolver & non_linear_solver,
                       SolverCallback & solver_callback) = 0;

  /// get instance of a time step solver
  virtual TimeStepSolver & getTimeStepSolver(const ID & time_step_solver_id);

  /// check if the given solver exists
  bool hasTimeStepSolver(const ID & solver_id) const;

  /* ------------------------------------------------------------------------ */
  const Mesh & getMesh() {
    if (mesh != nullptr) {
      return *mesh;
    }
    AKANTU_EXCEPTION("No mesh registered in this dof manager");
  }

  /* ------------------------------------------------------------------------ */
  AKANTU_GET_MACRO(Communicator, communicator, const auto &);
  AKANTU_GET_MACRO_NOT_CONST(Communicator, communicator, auto &);

  /* ------------------------------------------------------------------------ */
  AKANTU_GET_MACRO(Solution, *(solution.get()), const auto &);
  AKANTU_GET_MACRO_NOT_CONST(Solution, *(solution.get()), auto &);

  AKANTU_GET_MACRO(Residual, *(residual.get()), const auto &);
  AKANTU_GET_MACRO_NOT_CONST(Residual, *(residual.get()), auto &);

  /* ------------------------------------------------------------------------ */
  /* MeshEventHandler interface                                               */
  /* ------------------------------------------------------------------------ */
protected:
  friend class GlobalDOFInfoDataAccessor;
  /// helper function for the DOFManager::onNodesAdded method
  virtual std::pair<UInt, UInt> updateNodalDOFs(const ID & dof_id,
                                                const Array<UInt> & nodes_list);

  template <typename Func>
  auto countDOFsForNodes(const DOFData & dof_data, UInt nb_nodes,
                         Func && getNode);

  void updateDOFsData(DOFData & dof_data, UInt nb_new_local_dofs,
                      UInt nb_new_pure_local, UInt nb_nodes,
                      const std::function<UInt(UInt)> & getNode);

  void updateDOFsData(DOFData & dof_data, UInt nb_new_local_dofs,
                      UInt nb_new_pure_local);

  auto computeFirstDOFIDs(UInt nb_new_local_dofs, UInt nb_new_pure_local);

  /// resize all the global information and takes the needed measure like
  /// cleaning matrices profiles
  virtual void resizeGlobalArrays();

public:
  /// function to implement to react on  akantu::NewNodesEvent
  void onNodesAdded(const Array<UInt> & nodes_list,
                    const NewNodesEvent & event) override;
  /// function to implement to react on  akantu::RemovedNodesEvent
  void onNodesRemoved(const Array<UInt> & nodes_list,
                      const Array<UInt> & new_numbering,
                      const RemovedNodesEvent & event) override;
  /// function to implement to react on  akantu::NewElementsEvent
  void onElementsAdded(const Array<Element> & elements_list,
                       const NewElementsEvent & event) override;
  /// function to implement to react on  akantu::RemovedElementsEvent
  void onElementsRemoved(const Array<Element> & elements_list,
                         const ElementTypeMapArray<UInt> & new_numbering,
                         const RemovedElementsEvent & event) override;
  /// function to implement to react on  akantu::ChangedElementsEvent
  void onElementsChanged(const Array<Element> & old_elements_list,
                         const Array<Element> & new_elements_list,
                         const ElementTypeMapArray<UInt> & new_numbering,
                         const ChangedElementsEvent & event) override;

protected:
  inline DOFData & getDOFData(const ID & dof_id);
  inline const DOFData & getDOFData(const ID & dof_id) const;
  template <class DOFData_>
  inline DOFData_ & getDOFDataTyped(const ID & dof_id);
  template <class DOFData_>
  inline const DOFData_ & getDOFDataTyped(const ID & dof_id) const;

  virtual std::unique_ptr<DOFData> getNewDOFData(const ID & dof_id) = 0;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// dof representations in the dof manager
  struct DOFData {
    DOFData() = delete;
    explicit DOFData(const ID & dof_id);
    virtual ~DOFData();

    /// DOF support type (nodal, general) this is needed to determine how the
    /// dof are shared among processors
    DOFSupportType support_type;

    ID group_support;

    /// Degree of freedom array
    Array<Real> * dof{nullptr};

    /// Blocked degree of freedoms array
    Array<bool> * blocked_dofs{nullptr};

    /// Degree of freedoms increment
    Array<Real> * increment{nullptr};

    /// Degree of freedoms at previous step
    Array<Real> * previous{nullptr};

    /// Solution associated to the dof
    Array<Real> solution;

    /* ---------------------------------------------------------------------- */
    /* data for dynamic simulations                                           */
    /* ---------------------------------------------------------------------- */
    /// Degree of freedom derivatives arrays
    std::vector<Array<Real> *> dof_derivatives;

    /* ---------------------------------------------------------------------- */
    /// number of dofs to consider locally for this dof id
    UInt local_nb_dofs{0};

    /// Number of purely local dofs
    UInt pure_local_nb_dofs{0};

    /// number of ghost dofs
    UInt ghosts_nb_dofs{0};

    /// local numbering equation numbers
    Array<Int> local_equation_number;

    /// associated node for _dst_nodal dofs only
    Array<UInt> associated_nodes;

    virtual Array<Int> & getLocalEquationsNumbers() {
      return local_equation_number;
    }
  };

  /// type to store dofs information
  using DOFStorage = std::map<ID, std::unique_ptr<DOFData>>;

  /// type to store all the matrices
  using SparseMatricesMap = std::map<ID, std::unique_ptr<SparseMatrix>>;

  /// type to store all the lumped matrices
  using LumpedMatricesMap = std::map<ID, std::unique_ptr<SolverVector>>;

  /// type to store all the non linear solver
  using NonLinearSolversMap = std::map<ID, std::unique_ptr<NonLinearSolver>>;

  /// type to store all the time step solver
  using TimeStepSolversMap = std::map<ID, std::unique_ptr<TimeStepSolver>>;

  /// store a reference to the dof arrays
  DOFStorage dofs;

  /// list of sparse matrices that where created
  SparseMatricesMap matrices;

  /// list of lumped matrices
  LumpedMatricesMap lumped_matrices;

  /// non linear solvers storage
  NonLinearSolversMap non_linear_solvers;

  /// time step solvers storage
  TimeStepSolversMap time_step_solvers;

  /// reference to the underlying mesh
  Mesh * mesh{nullptr};

  /// Total number of degrees of freedom (size with the ghosts)
  UInt local_system_size{0};

  /// Number of purely local dofs (size without the ghosts)
  UInt pure_local_system_size{0};

  /// Total number of degrees of freedom
  UInt system_size{0};

  /// rhs to the system of equation corresponding to the residual linked to the
  /// different dofs
  std::unique_ptr<SolverVector> residual;

  /// solution of the system of equation corresponding to the different dofs
  std::unique_ptr<SolverVector> solution;

  /// a vector that helps internally to perform some tasks
  std::unique_ptr<SolverVector> data_cache;

  /// define the dofs type, local, shared, ghost
  Array<NodeFlag> dofs_flag;

  /// equation number in global numbering
  Array<Int> global_equation_number;

  using equation_numbers_map = std::unordered_map<Int, Int>;

  /// dual information of global_equation_number
  equation_numbers_map global_to_local_mapping;

  /// Communicator used for this manager, should be the same as in the mesh if a
  /// mesh is registered
  Communicator & communicator;

  /// accumulator to know what would be the next global id to use
  UInt first_global_dof_id{0};

  /// Release at last apply boundary on jacobian
  UInt jacobian_release{0};

  /// blocked degree of freedom in the system equation corresponding to the
  /// different dofs
  Array<Int> global_blocked_dofs;

  UInt global_blocked_dofs_release{0};

  /// blocked degree of freedom in the system equation corresponding to the
  /// different dofs
  Array<Int> previous_global_blocked_dofs;

  UInt previous_global_blocked_dofs_release{0};

private:
  /// This is for unit testing
  friend class DOFManagerTester;
};

using DefaultDOFManagerFactory =
    Factory<DOFManager, ID, const ID &, const MemoryID &>;
using DOFManagerFactory =
    Factory<DOFManager, ID, Mesh &, const ID &, const MemoryID &>;

} // namespace akantu

#include "dof_manager_inline_impl.hh"

#endif /* AKANTU_DOF_MANAGER_HH_ */
