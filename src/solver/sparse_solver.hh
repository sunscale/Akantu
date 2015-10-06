/**
 * @file   sparse_solver.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Aurelia Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 * @date creation: Mon Dec 13 2010
 * @date last modification: Mon Sep 15 2014
 *
 * @brief  interface for solvers
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef __AKANTU_SOLVER_HH__
#define __AKANTU_SOLVER_HH__

/* -------------------------------------------------------------------------- */
#include "static_solver.hh"
#include "data_accessor.hh"
#include "parsable.hh"

/* -------------------------------------------------------------------------- */

enum SolverParallelMethod {
  _not_parallel,
  _fully_distributed,
  _master_slave_distributed
};

namespace akantu {
  class DOFManager;
}

__BEGIN_AKANTU__

class SparseSolver : protected Memory,
                     public Parsable,
                     public StaticSolverEventHandler,
                     public DataAccessor {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  SparseSolver(DOFManager & dof_manager, const ID & matrix_id,
               const ID & id = "solver", const MemoryID & memory_id = 0);

  virtual ~SparseSolver();
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
};

__END_AKANTU__

#endif /* __AKANTU_SOLVER_HH__ */
