/**
 * @file   solver.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
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
#include "aka_common.hh"
#include "aka_memory.hh"
#include "aka_array.hh"
#include "sparse_matrix.hh"
#include "mesh.hh"
#include "static_communicator.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class SolverOptions {
public:
  SolverOptions(bool no_option = false) : no_option(no_option) { }

  virtual ~SolverOptions() {}

private:
  bool no_option;
};

extern SolverOptions _solver_no_options;

class Solver : protected Memory {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  Solver(SparseMatrix & matrix,
	 const ID & id = "solver",
	 const MemoryID & memory_id = 0);

  virtual ~Solver();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// initialize the solver
  virtual void initialize(SolverOptions & options = _solver_no_options) = 0;

  virtual void analysis() {};

  virtual void factorize() {};

  /// solve
  virtual void solve(Array<Real> & solution) = 0;
  virtual void solve() = 0;

  virtual void setRHS(const Array<Real> & rhs) = 0;

  /// function to print the contain of the class
  //  virtual void printself(std::ostream & stream, int indent = 0) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  AKANTU_GET_MACRO(RHS, *rhs, Array<Real> &);


  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:

  /// the matrix
  SparseMatrix * matrix;

  /// is sparse matrix allocated
  bool is_matrix_allocated;

  /// the right hand side
  Array<Real> * rhs;

  /// mesh
  const Mesh * mesh;

  /// pointer to the communicator
  StaticCommunicator & communicator;
};


__END_AKANTU__

#endif /* __AKANTU_SOLVER_HH__ */
