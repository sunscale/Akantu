/**
 * @file   static_solver.hh
 * @author Aurelia Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Wed Jul 30 14:44:12 2014
 *
 * @brief  Class handeling the initialization of external solvers
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
#include "aka_common.hh"


#ifndef __AKANTU_STATIC_SOLVER_HH__
#define __AKANTU_STATIC_SOLVER_HH__

__BEGIN_AKANTU__

class StaticSolver {
  /* ------------------------------------------------------------------------ */
  /* Constructors                                                             */
  /* ------------------------------------------------------------------------ */
private:
  StaticSolver() : is_initialized(false) {};

public:
  ~StaticSolver();

  /* ------------------------------------------------------------------------ */
  /// get an instance to the static solver
  static StaticSolver & getStaticSolver();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:  
  /// initialize what is needed for the compiled solver interfaces
  void initialize(int & argc, char ** & argv);

  /// finalize what is needed for the compiled solver interfaces
  void finalize();

  /* ------------------------------------------------------------------------ */
  /* Members                                                                  */
  /* ------------------------------------------------------------------------ */
private:
  bool is_initialized;
  static UInt nb_references;
  static StaticSolver * static_solver;
};

__END_AKANTU__


#endif /* __AKANTU_STATIC_SOLVER_HH__ */
