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
#include "aka_factory.hh"
#include "aka_memory.hh"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */
#include <map>
#include <set>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_DOF_MANAGER_HH__
#define __AKANTU_DOF_MANAGER_HH__

namespace akantu {
class TermsToAssemble;
class NonLinearSolver;
class TimeStepSolver;
class SparseMatrix;
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
private:
  /// common function to help registering dofs
  void registerDOFsInternal(const ID & dof_id, Array<Real> & dofs_array);

public:
  /// register an array of degree of freedom
  virtual void registerDOFs(const ID & dof_id, Array<Real> & dofs_array,
                            const DOFSupportType & support_type);

  /// the dof as an implied type of _dst_nodal and is defined only on a subset
  /// of nodes
  virtual void registerDOFs(const ID & dof_id, Array<Real> & dofs_array,
                            const ID & group_support);

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
                                  Real scale_factor = 1.) = 0;

  /// Assemble an array to the global lumped matrix array
  virtual void assembleToLumpedMatrix(const ID & dof_id,
                                      Array<Real> & array_to_assemble,
                                      const ID & lumped_mtx,
                                      Real scale_factor = 1.) = 0;

  /**
   * Assemble elementary values to a local array of the size nb_nodes *
   * nb_dof_per_node. The dof number is implicitly considered as
   * conn(el, n) * nb_nodes_per_element + d.
   * With 0 < n < nb_nodes_per_element and 0 < d < nb_dof_per_node
   **/
  virtual void assembleElementalArrayLocalArray(
      const Array<Real> & elementary_vect, Array<Real> & array_assembeled,
      const ElementType & type, const GhostType & ghost_type,
      Real scale_factor = 1.,
      const Array<UInt> & filter_elements = empty_filter);

  /**
   * Assemble elementary values to the global residual array. The dof number is
   * implicitly considered as conn(el, n) * nb_nodes_per_element + d.
   * With 0 < n < nb_nodes_per_element and 0 < d < nb_dof_per_node
   **/
  virtual void assembleElementalArrayToResidual(
      const ID & dof_id, const Array<Real> & elementary_vect,
      const ElementType & type, const GhostType & ghost_type,
      Real scale_factor = 1.,
      const Array<UInt> & filter_elements = empty_filter);

  /**
   * Assemble elementary values to a global array corresponding to a lumped
   * matrix
   */
  virtual void assembleElementalArrayToLumpedMatrix(
      const ID & dof_id, const Array<Real> & elementary_vect,
      const ID & lumped_mtx, const ElementType & type,
      const GhostType & ghost_type, Real scale_factor = 1.,
      const Array<UInt> & filter_elements = empty_filter);

  /**
   * Assemble elementary values to the global residual array. The dof number is
   * implicitly considered as conn(el, n) * nb_nodes_per_element + d.  With 0 <
   * n < nb_nodes_per_element and 0 < d < nb_dof_per_node
   **/
  virtual void assembleElementalMatricesToMatrix(
      const ID & matrix_id, const ID & dof_id,
      const Array<Real> & elementary_mat, const ElementType & type,
      const GhostType & ghost_type = _not_ghost,
      const MatrixType & elemental_matrix_type = _symmetric,
      const Array<UInt> & filter_elements = empty_filter) = 0;

  /// multiply a vector by a matrix and assemble the result to the residual
  virtual void assembleMatMulVectToResidual(const ID & dof_id, const ID & A_id,
                                            const Array<Real> & x,
                                            Real scale_factor = 1) = 0;

  /// multiply the dofs by a matrix and assemble the result to the residual
  virtual void assembleMatMulDOFsToResidual(const ID & A_id,
                                            Real scale_factor = 1);

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

  /// sets the residual to 0
  virtual void clearResidual() = 0;
  /// sets the matrix to 0
  virtual void clearMatrix(const ID & mtx) = 0;
  /// sets the lumped matrix to 0
  virtual void clearLumpedMatrix(const ID & mtx) = 0;

  /// splits the solution storage from a global view to the per dof storages
  void splitSolutionPerDOFs();

  /// extract a lumped matrix part corresponding to a given dof
  virtual void getLumpedMatrixPerDOFs(const ID & dof_id, const ID & lumped_mtx,
                                      Array<Real> & lumped) = 0;

protected:
  /// minimum functionality to implement per derived version of the DOFManager
  /// to allow the splitSolutionPerDOFs function to work
  virtual void getSolutionPerDOFs(const ID & dof_id,
                                  Array<Real> & solution_array) = 0;

protected:
  /* ------------------------------------------------------------------------ */
  /// register a matrix
  SparseMatrix & registerSparseMatrix(const ID & matrix_id,
                                      std::unique_ptr<SparseMatrix> & matrix);

  /// register a non linear solver instantiated by a derived class
  NonLinearSolver &
  registerNonLinearSolver(const ID & non_linear_solver_id,
                          std::unique_ptr<NonLinearSolver> & non_linear_solver);

  /// register a time step solver instantiated by a derived class
  TimeStepSolver &
  registerTimeStepSolver(const ID & time_step_solver_id,
                         std::unique_ptr<TimeStepSolver> & time_step_solver);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// Get the equation numbers corresponding to a dof_id. This might be used to
  /// access the matrix.
  inline const Array<UInt> & getEquationsNumbers(const ID & dof_id) const;

  /// Global number of dofs
  AKANTU_GET_MACRO(SystemSize, this->system_size, UInt);

  /// Local number of dofs
  AKANTU_GET_MACRO(LocalSystemSize, this->local_system_size, UInt);

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
  inline bool hasDOFs(const ID & dofs_id) const;

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

  /// Get the reference of an existing matrix
  SparseMatrix & getMatrix(const ID & matrix_id);

  /// check if the given matrix exists
  bool hasMatrix(const ID & matrix_id) const;

  /// Get an instance of a new lumped matrix
  virtual Array<Real> & getNewLumpedMatrix(const ID & matrix_id);
  /// Get the lumped version of a given matrix
  const Array<Real> & getLumpedMatrix(const ID & matrix_id) const;
  /// Get the lumped version of a given matrix
  Array<Real> & getLumpedMatrix(const ID & matrix_id);

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
                       NonLinearSolver & non_linear_solver) = 0;

  /// get instance of a time step solver
  virtual TimeStepSolver & getTimeStepSolver(const ID & time_step_solver_id);

  /// check if the given solver exists
  bool hasTimeStepSolver(const ID & solver_id) const;

  /* ------------------------------------------------------------------------ */
  const Mesh & getMesh() {
    if (mesh) {
      return *mesh;
    } else {
      AKANTU_EXCEPTION("No mesh registered in this dof manager");
    }
  }

  /* ------------------------------------------------------------------------ */
  AKANTU_GET_MACRO(Communicator, communicator, const auto &);
  AKANTU_GET_MACRO_NOT_CONST(Communicator, communicator, auto &);

  /* ------------------------------------------------------------------------ */
  /* MeshEventHandler interface                                               */
  /* ------------------------------------------------------------------------ */
protected:
  /// helper function for the DOFManager::onNodesAdded method
  virtual std::pair<UInt, UInt> updateNodalDOFs(const ID & dof_id,
                                                const Array<UInt> & nodes_list);

  template <typename Func>
  auto countDOFsForNodes(const DOFData & dof_data, UInt nb_nodes,
                        Func && getNode);

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
  template <class _DOFData>
  inline _DOFData & getDOFDataTyped(const ID & dof_id);
  template <class _DOFData>
  inline const _DOFData & getDOFDataTyped(const ID & dof_id) const;

  virtual DOFData & getNewDOFData(const ID & dof_id);

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
    Array<Real> * dof;

    /// Blocked degree of freedoms array
    Array<bool> * blocked_dofs;

    /// Degree of freedoms increment
    Array<Real> * increment;

    /// Degree of freedoms at previous step
    Array<Real> * previous;

    /// Solution associated to the dof
    Array<Real> solution;

    /// local numbering equation numbers
    Array<UInt> local_equation_number;

    /* ---------------------------------------------------------------------- */
    /* data for dynamic simulations                                           */
    /* ---------------------------------------------------------------------- */
    /// Degree of freedom derivatives arrays
    std::vector<Array<Real> *> dof_derivatives;
  };

  /// type to store dofs information
  using DOFStorage = std::map<ID, std::unique_ptr<DOFData>>;

  /// type to store all the matrices
  using SparseMatricesMap = std::map<ID, std::unique_ptr<SparseMatrix>>;

  /// type to store all the lumped matrices
  using LumpedMatricesMap = std::map<ID, std::unique_ptr<Array<Real>>>;

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

  /// Communicator used for this manager, should be the same as in the mesh if a
  /// mesh is registered
  Communicator & communicator;
};

using DefaultDOFManagerFactory =
    Factory<DOFManager, ID, const ID &, const MemoryID &>;
using DOFManagerFactory =
    Factory<DOFManager, ID, Mesh &, const ID &, const MemoryID &>;

} // namespace akantu

#include "dof_manager_inline_impl.cc"

#endif /* __AKANTU_DOF_MANAGER_HH__ */
