/**
 * @file   solver_petsc.hh
 *
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue May 13 2014
 * @date last modification: Mon Jun 19 2017
 *
 * @brief  Solver class interface for the petsc solver
 *
 *
 * Copyright (©) 2014-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "sparse_solver.hh"
/* -------------------------------------------------------------------------- */
#include <petscksp.h>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_SOLVER_PETSC_HH__
#define __AKANTU_SOLVER_PETSC_HH__

namespace akantu {
class SparseMatrixPETSc;
class DOFManagerPETSc;
} // namespace akantu

namespace akantu {

class SolverPETSc : public SparseSolver {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  SolverPETSc(DOFManagerPETSc & dof_manager, const ID & matrix_id,
              const ID & id = "solver_petsc", const MemoryID & memory_id = 0);

  virtual ~SolverPETSc();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// create the solver context and set the matrices
  virtual void setOperators();
  virtual void solve();

private:
  /// DOFManager correctly typed
  DOFManagerPETSc & dof_manager;

  /// PETSc linear solver
  KSP ksp;

  /// Matrix defining the system of equations
  SparseMatrixPETSc & matrix;

  /// specify if the petsc_data is initialized or not
  bool is_petsc_data_initialized;
};

//   SolverPETSc(int argc, char *argv[]) : allocated_(false) {

//     /*
//      Set linear solver defaults for this problem (optional).
//      - By extracting the KSP and PC contexts from the KSP context,
//      we can then directly call any KSP and PC routines to set
//      various options.
//      - The following four statements are optional; all of these
//      parameters could alternatively be specified at runtime via
//      KSPSetFromOptions();
//      */
//     //      ierr = KSPGetPC(ksp_,&pc);CHKERRCONTINUE(ierr);
//     //      ierr = PCSetType(pc,PCILU);CHKERRCONTINUE(ierr);
//     //    ierr = PCSetType(pc,PCJACOBI);CHKERRCONTINUE(ierr);
//     ierr =
//     KSPSetTolerances(ksp_,1.e-5,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRCONTINUE(ierr);
//   }

//   //! Overload operator() to solve system of linear equations
//   sparse_vector_type operator()(const sparse_matrix_type& AA, const
//   sparse_vector_type& bb);

//   //! Overload operator() to obtain reaction vector
//   sparse_vector_type operator()(const sparse_matrix_type& Kpf, const
//   sparse_matrix_type& Kpp, const sparse_vector_type& Up);

//   //! Overload operator() to obtain the addition two vectors
//   sparse_vector_type operator()(const sparse_vector_type& aa, const
//   sparse_vector_type& bb);

//   value_type norm(const sparse_matrix_type& aa, Element_insertion_type it =
//   Add_t);

//   value_type norm(const sparse_vector_type& aa, Element_insertion_type it =
//   Add_t);

//   // NOTE: the destructor will return an error if it is called after
//   MPI_Finalize is
//   // called because it uses collect communication to free-up allocated
//   memory.
//   ~SolverPETSc() {

//     static bool exit = false;
//     if (!exit) {
//       // add finalize PETSc function at exit
//       atexit(finalize);
//       exit = true;
//     }

//     if (allocated_) {
//       PetscErrorCode ierr = MatDestroy(&A_);CHKERRCONTINUE(ierr);
//       ierr = VecDestroy(&x_);CHKERRCONTINUE(ierr);
//       ierr = KSPDestroy(&ksp_);CHKERRCONTINUE(ierr);
//     }
//   }

//   /* from the PETSc library, these are the options that can be passed
//    to the command line

//    Options Database Keys

//    -options_table	                - Calls PetscOptionsView()
//    -options_left	                - Prints unused options that remain in
//    the
//    database
//    -objects_left                  - Prints list of all objects that have not
//    been freed
//    -mpidump	                    - Calls PetscMPIDump()
//    -malloc_dump	                - Calls PetscMallocDump()
//    -malloc_info	                - Prints total memory usage
//    -malloc_log	                - Prints summary of memory usage

//    Options Database Keys for Profiling

//    -log_summary [filename]	    - Prints summary of flop and timing
//    information to screen.
//    If the filename is specified the summary is written to the file. See
//    PetscLogView().
//    -log_summary_python [filename]	- Prints data on of flop and timing
//    usage
//    to a file or screen.
//    -log_all [filename]	        - Logs extensive profiling information
//    See
//    PetscLogDump().
//    -log [filename]	            - Logs basic profiline information See
//    PetscLogDump().
//    -log_sync	                    - Log the synchronization in scatters,
//    inner products and norms
//    -log_mpe [filename]            - Creates a logfile viewable by the utility
//    Upshot/Nupshot (in MPICH distribution)
//     }
//   }
// };

} // namespace akantu

#endif /* __AKANTU_SOLVER_PETSC_HH__ */
