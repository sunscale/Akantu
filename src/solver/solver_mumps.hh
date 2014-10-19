/**
 * @file   solver_mumps.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Dec 13 2010
 * @date last modification: Mon Sep 15 2014
 *
 * @brief  Solver class implementation for the mumps solver
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
#ifndef __AKANTU_SOLVER_MUMPS_HH__
#define __AKANTU_SOLVER_MUMPS_HH__
#include <dmumps_c.h>

#include "solver.hh"
#include "static_communicator.hh"

__BEGIN_AKANTU__

class SolverMumpsOptions : public SolverOptions {
public:
  enum ParallelMethod {
    _not_parallel,
    _fully_distributed,
    _master_slave_distributed
  };

  SolverMumpsOptions(ParallelMethod parallel_method = _fully_distributed) :
    SolverOptions(),
    parallel_method(parallel_method) { }

private:
  friend class SolverMumps;
  ParallelMethod parallel_method;
};

class SolverMumps : public Solver, public CommunicatorEventHandler {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  SolverMumps(SparseMatrix & sparse_matrix,
	      const ID & id = "solver_mumps",
	      const MemoryID & memory_id = 0);

  virtual ~SolverMumps();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// build the profile and do the analysis part
  void initialize(SolverOptions & options = _solver_no_options);

  void initializeSlave(SolverOptions & options = _solver_no_options);

  /// analysis (symbolic facto + permutations)
  void analysis();

  /// factorize the matrix
  void factorize();

  /// solve the system
  void solve(Array<Real> & solution);
  void solve();

  void solveSlave();

  virtual void setRHS(const Array<Real> & rhs);


  virtual void onCommunicatorFinalize(const StaticCommunicator & communicator);

private:
  /// print the error if any happened in mumps
  void printError();

  /// clean the mumps info
  void destroyMumpsData();

  /// set internal values;
  void initMumpsData();

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
private:
  /// access the control variable
  inline Int & icntl(UInt i) {
    return mumps_data.icntl[i - 1];
  }

  /// access the results info
  inline Int & info(UInt i) {
    return mumps_data.info[i - 1];
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  /// mumps data
  DMUMPS_STRUC_C mumps_data;

  /// specify if the mumps_data are initialized or not
  bool is_mumps_data_initialized;

  UInt prank;

  /* ------------------------------------------------------------------------ */
  /* Local types                                                              */
  /* ------------------------------------------------------------------------ */
private:
  SolverMumpsOptions::ParallelMethod parallel_method;

  bool rhs_is_local;

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
