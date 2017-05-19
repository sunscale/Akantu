/**
 * @file   dof_manager.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Wed Jul 22 11:43:43 2015
 *
 * @brief  Class handling the different types of dofs
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
}

namespace akantu {

class DOFManager : protected Memory, protected MeshEventHandler {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  DOFManager(const ID & id = "dof_manager", const MemoryID & memory_id = 0);
  DOFManager(Mesh & mesh, const ID & id = "dof_manager", const MemoryID & memory_id = 0);
  virtual ~DOFManager();

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
                                  const Array<Real> & array_to_assemble,
                                  Real scale_factor = 1.) = 0;

  /// Assemble an array to the global lumped matrix array
  virtual void assembleToLumpedMatrix(const ID & dof_id,
                                      const Array<Real> & array_to_assemble,
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
  void registerSparseMatrix(const ID & matrix_id, SparseMatrix & matrix);

  /// register a non linear solver instantiated by a derived class
  void registerNonLinearSolver(const ID & non_linear_solver_id,
                               NonLinearSolver & non_linear_solver);

  /// register a time step solver instantiated by a derived class
  void registerTimeStepSolver(const ID & time_step_solver_id,
                              TimeStepSolver & time_step_solver);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// Get the equation numbers corresponding to a dof_id. This might be used to
  /// access the matrix.
  virtual void getEquationsNumbers(const ID & dof_id,
                                   Array<UInt> & equation_numbers);

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
  AKANTU_GET_MACRO(Communicator, communicator, const StaticCommunicator &);

  /* ------------------------------------------------------------------------ */
  /* MeshEventHandler interface                                               */
  /* ------------------------------------------------------------------------ */
private:
  /// fills the nodes_to_elements structure
  void fillNodesToElements();

public:
  /// function to implement to react on  akantu::NewNodesEvent
  virtual void onNodesAdded(const Array<UInt> & nodes_list,
                            const NewNodesEvent & event);
  /// function to implement to react on  akantu::RemovedNodesEvent
  virtual void onNodesRemoved(const Array<UInt> & nodes_list,
                              const Array<UInt> & new_numbering,
                              const RemovedNodesEvent & event);
  /// function to implement to react on  akantu::NewElementsEvent
  virtual void onElementsAdded(const Array<Element> & elements_list,
                               const NewElementsEvent & event);
  /// function to implement to react on  akantu::RemovedElementsEvent
  virtual void
  onElementsRemoved(const Array<Element> & elements_list,
                    const ElementTypeMapArray<UInt> & new_numbering,
                    const RemovedElementsEvent & event);
  /// function to implement to react on  akantu::ChangedElementsEvent
  virtual void
  onElementsChanged(const Array<Element> & old_elements_list,
                    const Array<Element> & new_elements_list,
                    const ElementTypeMapArray<UInt> & new_numbering,
                    const ChangedElementsEvent & event);

protected:
  struct DOFData;
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

  typedef Array<std::set<Element> *> NodesToElements;

  /// This info is stored to simplify the dynamic changes
  NodesToElements nodes_to_elements;

  /// type to store dofs information
  typedef std::map<ID, DOFData *> DOFStorage;

  /// type to store all the matrices
  typedef std::map<ID, SparseMatrix *> SparseMatricesMap;

  /// type to store all the lumped matrices
  typedef std::map<ID, Array<Real> *> LumpedMatricesMap;

  /// type to store all the non linear solver
  typedef std::map<ID, NonLinearSolver *> NonLinearSolversMap;

  /// type to store all the time step solver
  typedef std::map<ID, TimeStepSolver *> TimeStepSolversMap;

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
  Mesh * mesh;

  /// Total number of degrees of freedom (size with the ghosts)
  UInt local_system_size;

  /// Number of purely local dofs (size without the ghosts)
  UInt pure_local_system_size;

  /// Total number of degrees of freedom
  UInt system_size;

  /// Communicator used for this manager, should be the same as in the mesh if a
  /// mesh is registered
  const StaticCommunicator & communicator;
};

} // akantu

#include "dof_manager_inline_impl.cc"

#endif /* __AKANTU_DOF_MANAGER_HH__ */
