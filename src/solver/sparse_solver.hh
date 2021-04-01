/**
 * @file   sparse_solver.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Wed Jan 24 2018
 *
 * @brief  interface for solvers
 *
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "communicator_event_handler.hh"
#include "parsable.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_SOLVER_HH_
#define AKANTU_SOLVER_HH_

namespace akantu {
enum SolverParallelMethod {
  _not_parallel,
  _fully_distributed,
  _master_slave_distributed
};

class DOFManager;
} // namespace akantu

namespace akantu {

class SparseSolver : public Parsable, public CommunicatorEventHandler {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  SparseSolver(DOFManager & dof_manager, const ID & matrix_id,
               const ID & id = "solver");

  ~SparseSolver() override;
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// initialize the solver
  virtual void initialize() = 0;

  virtual void analysis(){};

  virtual void factorize(){};

  virtual void solve(){};

protected:
  virtual void destroyInternalData(){};

public:
  virtual void beforeStaticSolverDestroy();

  void createSynchronizerRegistry();
  /* ------------------------------------------------------------------------ */
  /* Data Accessor inherited members                                          */
  /* ------------------------------------------------------------------------ */
public:
  void onCommunicatorFinalize() override;

  // inline virtual UInt getNbDataForDOFs(const Array<UInt> & dofs,
  //                                      SynchronizationTag tag) const;

  // inline virtual void packDOFData(CommunicationBuffer & buffer,
  //                                 const Array<UInt> & dofs,
  //                                 SynchronizationTag tag) const;

  // inline virtual void unpackDOFData(CommunicationBuffer & buffer,
  //                                   const Array<UInt> & dofs,
  //                                   SynchronizationTag tag);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// manager handling the dofs for this SparseMatrix solver
  DOFManager & _dof_manager;

  /// The id of the associated matrix
  ID matrix_id;

  /// How to parallelize the solve
  SolverParallelMethod parallel_method;

  /// Communicator used by the solver
  Communicator & communicator;
};

namespace debug {
  class SingularMatrixException : public Exception {
  public:
    SingularMatrixException(const SparseMatrix & matrix)
        : Exception("Solver encountered singular matrix"), matrix(matrix) {}
    const SparseMatrix & matrix;
  };
} // namespace debug

} // namespace akantu

#endif /* AKANTU_SOLVER_HH_ */
