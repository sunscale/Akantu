/**
 * @file   sparse_solver_mumps.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Tue Jan 19 2016
 *
 * @brief  Solver class implementation for the mumps solver
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include "sparse_solver.hh"
/* -------------------------------------------------------------------------- */
#include <dmumps_c.h>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_SOLVER_MUMPS_HH__
#define __AKANTU_SOLVER_MUMPS_HH__

namespace akantu {
  class DOFManagerDefault;
  class SparseMatrixAIJ;
}


__BEGIN_AKANTU__

class SparseSolverMumps : public SparseSolver {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  SparseSolverMumps(DOFManagerDefault & dof_manager,
                    const ID & matrix_id,
                    const ID & id = "sparse_solver_mumps",
                    const MemoryID & memory_id = 0);

  virtual ~SparseSolverMumps();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// build the profile and do the analysis part
  virtual void initialize();

  //  void initializeSlave(SolverOptions & options = _solver_no_options);

  /// analysis (symbolic facto + permutations)
  virtual void analysis();

  /// factorize the matrix
  virtual void factorize();

  /// solve the system
  virtual void solve();

private:
  /// print the error if any happened in mumps
  void printError();

  /// set internal values;
  void initMumpsData();

  /// set the level of verbosity of mumps based on the debug level of akantu
  void setOutputLevel();

protected:
  /// de-initialize the internal data
  virtual void destroyInternalData();

  /// check if initialized and except if it is not the case
  void checkInitialized();

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
private:
  /// access the control variable
  inline Int & icntl(UInt i) { return mumps_data.icntl[i - 1]; }

  /// access the results info
  inline Int & info(UInt i) { return mumps_data.info[i - 1]; }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// DOFManager used by the Mumps implementation of the SparseSolver
  DOFManagerDefault & dof_manager;

  /// AIJ Matrix, usualy the jacobian matrix
  SparseMatrixAIJ & matrix;

  /// Right hand side per processors
  const Array<Real> & rhs;

  /// solution vector per processors
  Array<Real> & solution;

  /// Full right hand side on the master processors and solution after solve
  Array<Real> master_rhs_solution;

  /// mumps data
  DMUMPS_STRUC_C mumps_data;

  /// Rank of the current process
  UInt prank;

  /// matrix release at last solve
  UInt last_profile_release;

  /// matrix release at last solve
  UInt last_value_release;

  /// check if the solver data are initialized
  bool is_initialized;

  /* ------------------------------------------------------------------------ */
  /* Local types                                                              */
  /* ------------------------------------------------------------------------ */
private:
  SolverParallelMethod parallel_method;

  // bool rhs_is_local;

  enum SolverMumpsJob {
    _smj_initialize = -1,
    _smj_analyze = 1,
    _smj_factorize = 2,
    _smj_solve = 3,
    _smj_analyze_factorize = 4,
    _smj_factorize_solve = 5,
    _smj_complete = 6, // analyze, factorize, solve
    _smj_destroy = -2
  };
};

__END_AKANTU__

#endif /* __AKANTU_SOLVER_MUMPS_HH__ */
