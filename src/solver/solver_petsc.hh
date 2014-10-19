/**
 * @file   solver_petsc.hh
 *
 # @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Mon Dec 13 10:48:06 2010
 *
 * @brief  Solver class implementation for the petsc solver
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
// #ifndef __AKANTU_SOLVER_PETSC_HH__
// #define __AKANTU_SOLVER_PETSC_HH__


// #include <petscksp.h>


// #include "solver.hh"
// #include "static_communicator.hh"
// #include "sparse_matrix.hh"

// __BEGIN_AKANTU__



// struct SolverPETSc : public Solver, public CommunicatorEventHandler {
  
//   typedef double value_type;
//   typedef sparse_vector<value_type> sparse_vector_type;
//   typedef SparseMatrix sparse_matrix_type;
  
//   Mat A_;                //!< linear system matrix
//   Vec x_;                //!< Solution vector
//   KSP ksp_;              //!< linear solver context
  
//   bool allocated_;
  
//   SolverPETSc(int argc, char *argv[]) : allocated_(false) {
    
    
//     PetscInitialize(&argc, &argv,NULL,NULL);
//     PetscErrorCode ierr;
    
//     // create linear solver context
//     ierr = KSPCreate(PETSC_COMM_WORLD, &ksp_);CHKERRCONTINUE(ierr);
    
//     // initial nonzero guess
//     ierr = KSPSetInitialGuessNonzero(ksp_,PETSC_TRUE); CHKERRCONTINUE(ierr);
    
//     // set runtime options
//     ierr = KSPSetFromOptions(ksp_);CHKERRCONTINUE(ierr);
    
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
//     ierr = KSPSetTolerances(ksp_,1.e-5,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRCONTINUE(ierr);
//   }
  
//   //! Overload operator() to solve system of linear equations
//   sparse_vector_type operator()(const sparse_matrix_type& AA, const sparse_vector_type& bb);
  
//   //! Overload operator() to obtain reaction vector
//   sparse_vector_type operator()(const sparse_matrix_type& Kpf, const sparse_matrix_type& Kpp, const sparse_vector_type& Up);
  
//   //! Overload operator() to obtain the addition two vectors
//   sparse_vector_type operator()(const sparse_vector_type& aa, const sparse_vector_type& bb);
  
//   value_type norm(const sparse_matrix_type& aa, Element_insertion_type it = Add_t);
  
//   value_type norm(const sparse_vector_type& aa, Element_insertion_type it = Add_t);
  
//   // NOTE: the destructor will return an error if it is called after MPI_Finalize is
//   // called because it uses collect communication to free-up allocated memory.
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
//    -options_left	                - Prints unused options that remain in the database
//    -objects_left                  - Prints list of all objects that have not been freed
//    -mpidump	                    - Calls PetscMPIDump()
//    -malloc_dump	                - Calls PetscMallocDump()
//    -malloc_info	                - Prints total memory usage
//    -malloc_log	                - Prints summary of memory usage
   
//    Options Database Keys for Profiling
   
//    -log_summary [filename]	    - Prints summary of flop and timing information to screen.
//    If the filename is specified the summary is written to the file. See PetscLogView().
//    -log_summary_python [filename]	- Prints data on of flop and timing usage to a file or screen.
//    -log_all [filename]	        - Logs extensive profiling information See PetscLogDump().
//    -log [filename]	            - Logs basic profiline information See PetscLogDump().
//    -log_sync	                    - Log the synchronization in scatters, inner products and norms
//    -log_mpe [filename]            - Creates a logfile viewable by the utility Upshot/Nupshot (in MPICH distribution)
//    */
//   static void finalize() {
    
//     static bool finalized = false;
//     if (!finalized) {
      
//       cout<<"*** INFO *** PETSc is finalizing..."<<endl;
      
//       // finalize PETSc
//       PetscErrorCode ierr = PetscFinalize();CHKERRCONTINUE(ierr);
//       finalized = true;
      
//       cout<<"*** INFO *** Process "<<Parallel_base::rank_<<" is finalizing..."<<endl;
      
//       // finalize MPI
//       MPI_Finalize();
//     }
//   }
// };




// class SolverPETSc : public Solver, public CommunicatorEventHandler {
//   /* ------------------------------------------------------------------------ */
//   /* Constructors/Destructors                                                 */
//   /* ------------------------------------------------------------------------ */
// public:
  
//   SolverPETSc(SparseMatrix & sparse_matrix,
//               const ID & id = "solver_petsc",
//               const MemoryID & memory_id = 0);
  
//   virtual SolverPETSc();
  
//   /* ------------------------------------------------------------------------ */
//   /* Methods                                                                  */
//   /* ------------------------------------------------------------------------ */
// public:
  
//   /// build the profile and do the analysis part
//   void initialize(SolverOptions & options = _solver_no_options);
  
//   void initializeSlave(SolverOptions & options = _solver_no_options);
  
  
//   /// factorize and solve the system
//   void solve(Array<Real> & solution);
//   void solve();
  
//   void solveSlave();
  
//   virtual void setRHS(Array<Real> & rhs);
  
//   /// function to print the contain of the class
//   //  virtual void printself(std::ostream & stream, int indent = 0) const;
  
//   virtual void onCommunicatorFinalize(const StaticCommunicator & communicator);
  
// private:
  
//   void destroyMumpsData();
  
//   inline Int & icntl(UInt i) {
//     return mumps_data.icntl[i - 1];
//   }
  
//   inline Int & info(UInt i) {
//     return mumps_data.info[i - 1];
//   }
  
//   void initMumpsData(SolverMumpsOptions::ParallelMethod parallel_method);
  
//   /* ------------------------------------------------------------------------ */
//   /* Accessors                                                                */
//   /* ------------------------------------------------------------------------ */
// public:
  
//   /* ------------------------------------------------------------------------ */
//   /* Class Members                                                            */
//   /* ------------------------------------------------------------------------ */
// private:
  
//   /// mumps data
//   DMUMPS_STRUC_C mumps_data;
  
//   /// specify if the mumps_data are initialized or not
//   bool is_mumps_data_initialized;
  
//   UInt prank;
  
//   /* ------------------------------------------------------------------------ */
//   /* Local types                                                              */
//   /* ------------------------------------------------------------------------ */

// private:
//   SolverMumpsOptions::ParallelMethod parallel_method;
  
//   bool rhs_is_local;
  
//   enum SolverMumpsJob {
//     _smj_initialize = -1,
//     _smj_analyze = 1,
//     _smj_factorize = 2,
//     _smj_solve = 3,
//     _smj_analyze_factorize = 4,
//     _smj_factorize_solve = 5,
//     _smj_complete = 6, // analyze, factorize, solve
//     _smj_destroy = -2
//   };
// };


// /* -------------------------------------------------------------------------- */
// /* inline functions                                                           */
// /* -------------------------------------------------------------------------- */

// //#include "solver_mumps_inline_impl.cc"

// /// standard output stream operator
// // inline std::ostream & operator <<(std::ostream & stream, const SolverMumps & _this)
// // {
// //   _this.printself(stream);
// //   return stream;
// // }


// __END_AKANTU__

// #endif /* __AKANTU_SOLVER_PETSC_HH__ */
