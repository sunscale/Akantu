/**
 * @file   solver_petsc.cc
 *
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue May 13 2014
 * @date last modification: Sun Aug 13 2017
 *
 * @brief  Solver class implementation for the petsc solver
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
#include "solver_petsc.hh"
#include "dof_manager_petsc.hh"
#include "mpi_communicator_data.hh"
#include "solver_vector_petsc.hh"
#include "sparse_matrix_petsc.hh"
/* -------------------------------------------------------------------------- */
#include <petscksp.h>
//#include <petscsys.h>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
SolverPETSc::SolverPETSc(DOFManagerPETSc & dof_manager, const ID & matrix_id,
                         const ID & id, const MemoryID & memory_id)
    : SparseSolver(dof_manager, matrix_id, id, memory_id),
      dof_manager(dof_manager), matrix(dof_manager.getMatrix(matrix_id)),
      is_petsc_data_initialized(false) {
  auto mpi_comm = dof_manager.getMPIComm();

  /// create a solver context
  PETSc_call(KSPCreate, mpi_comm, &this->ksp);
}

/* -------------------------------------------------------------------------- */
SolverPETSc::~SolverPETSc() {
  AKANTU_DEBUG_IN();

  if (ksp)
    PETSc_call(KSPDestroy, &ksp);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolverPETSc::setOperators() {
  // set the matrix that defines the linear system and the matrix for
// preconditioning (here they are the same)
#if PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 5
  PETSc_call(KSPSetOperators, ksp, this->matrix.getMat(),
             this->matrix.getMat());
#else
  PETSc_call(KSPSetOperators, ksp, this->matrix.getMat(), this->matrix.getMat(),
             SAME_NONZERO_PATTERN);
#endif

  // If this is not called the solution vector is zeroed in the call to
  // KSPSolve().
  PETSc_call(KSPSetInitialGuessNonzero, ksp, PETSC_TRUE);
  PETSc_call(KSPSetFromOptions, ksp);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolverPETSc::solve() {
  Vec & rhs(this->dof_manager.getResidual());
  Vec & solution(this->dof_manager.getSolution());

  PETSc_call(KSPSolve, ksp, rhs, solution);
}

// /* --------------------------------------------------------------------------
// */
// void SolverPETSc::solve(Array<Real> & solution) {
//   AKANTU_DEBUG_IN();

//   this->solution = &solution;
//   this->solve();

//   PetscErrorCode ierr;
//   PETScMatrix * petsc_matrix = static_cast<PETScMatrix *>(this->matrix);

//   // ierr = VecGetOwnershipRange(this->petsc_solver_wrapper->solution,
//   // solution_begin, solution_end); CHKERRXX(ierr);
//   // ierr = VecGetArray(this->petsc_solver_wrapper->solution, PetscScalar
//   // **array); CHKERRXX(ierr);
//   UInt nb_component = solution.getNbComponent();
//   UInt size = solution.size();

//   for (UInt i = 0; i < size; ++i) {
//     for (UInt j = 0; j < nb_component; ++j) {
//       UInt i_local = i * nb_component + j;
//       if (this->matrix->getDOFSynchronizer().isLocalOrMasterDOF(i_local)) {
//         Int i_global =
//             this->matrix->getDOFSynchronizer().getDOFGlobalID(i_local);
//         ierr =
//         AOApplicationToPetsc(petsc_matrix->getPETScMatrixWrapper()->ao,
//                                     1, &(i_global));
//         CHKERRXX(ierr);
//         ierr = VecGetValues(this->petsc_solver_wrapper->solution, 1,
//         &i_global,
//                             &solution(i, j));
//         CHKERRXX(ierr);
//       }
//     }
//   }
//   synch_registry->synchronize(SynchronizationTag::_solver_solution);

//   AKANTU_DEBUG_OUT();
// }

/* -------------------------------------------------------------------------- */
// void finalize_petsc() {

//   static bool finalized = false;
//   if (!finalized) {

//     cout<<"*** INFO *** PETSc is finalizing..."<<endl;
//     // finalize PETSc
//     PetscErrorCode ierr = PetscFinalize();CHKERRCONTINUE(ierr);
//     finalized = true;
//   }
// }

// SolverPETSc::sparse_vector_type
// SolverPETSc::operator()(const SolverPETSc::sparse_matrix_type& AA,
//                          const SolverPETSc::sparse_vector_type& bb) {

// #ifdef CPPUTILS_VERBOSE
//   // parallel output stream
//   Output_stream out;
//   // timer
//   cpputils::ctimer timer;
//   out<<"Inside PETSc solver: "<<timer<<endl;
// #endif

// #ifdef CPPUTILS_VERBOSE
//   out<<"  Inside operator()(const sparse_matrix_type&, const
//   sparse_vector_type&)... "<<timer<<endl;
// #endif

//   assert(AA.rows() == bb.size());

//   //  KSP ksp;            //!< linear solver context

//   Vec            b;                /* RHS */
//   PC             pc;               /* preconditioner context */
//   PetscErrorCode ierr;
//   PetscInt       nlocal;
//   PetscInt       n = bb.size();
//   VecScatter ctx;

//   /*
//    Create vectors.  Note that we form 1 vector from scratch and
//    then duplicate as needed. For this simple case let PETSc decide how
//    many elements of the vector are stored on each processor. The second
//    argument to VecSetSizes() below causes PETSc to decide.
//    */
//   ierr = VecCreate(PETSC_COMM_WORLD,&b);CHKERRCONTINUE(ierr);
//   ierr = VecSetSizes(b,PETSC_DECIDE,n);CHKERRCONTINUE(ierr);
//   ierr = VecSetFromOptions(b);CHKERRCONTINUE(ierr);
//   if (!allocated_) {
//     ierr = VecDuplicate(b,&x_);CHKERRCONTINUE(ierr);
//   } else
//     VecZeroEntries(x_);

// #ifdef CPPUTILS_VERBOSE
//   out<<"  Vectors created... "<<timer<<endl;
// #endif

//   /* Set hight-hand-side vector */
//   for (sparse_vector_type::const_hash_iterator it = bb.map_.begin(); it !=
//   bb.map_.end(); ++it) {
//     int row = it->first;
//     ierr = VecSetValues(b, 1, &row, &it->second, ADD_VALUES);
//   }

// #ifdef CPPUTILS_VERBOSE
//   out<<"  Right hand side set... "<<timer<<endl;
// #endif

//   /*
//    Assemble vector, using the 2-step process:
//    VecAssemblyBegin(), VecAssemblyEnd()
//    Computations can be done while messages are in transition
//    by placing code between these two statements.
//    */
//   ierr = VecAssemblyBegin(b);CHKERRCONTINUE(ierr);
//   ierr = VecAssemblyEnd(b);CHKERRCONTINUE(ierr);

// #ifdef CPPUTILS_VERBOSE
//   out<<"  Right-hand-side vector assembled... "<<timer<<endl;

//   ierr = VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRCONTINUE(ierr);

//   Vec b_all;
//   ierr = VecScatterCreateToAll(b, &ctx, &b_all);CHKERRCONTINUE(ierr);
//   ierr =
//   VecScatterBegin(ctx,b,b_all,INSERT_VALUES,SCATTER_FORWARD);CHKERRCONTINUE(ierr);
//   ierr =
//   VecScatterEnd(ctx,b,b_all,INSERT_VALUES,SCATTER_FORWARD);CHKERRCONTINUE(ierr);

//   value_type nrm;
//   VecNorm(b_all,NORM_2,&nrm);
//   out<<"  Norm of right hand side... "<<nrm<<endl;
// #endif

//   /* Identify the starting and ending mesh points on each
//    processor for the interior part of the mesh. We let PETSc decide
//    above. */

//   //    PetscInt rstart,rend;
//   //    ierr = VecGetOwnershipRange(x_,&rstart,&rend);CHKERRCONTINUE(ierr);
//   ierr = VecGetLocalSize(x_,&nlocal);CHKERRCONTINUE(ierr);

//   /*
//    Create matrix.  When using MatCreate(), the matrix format can
//    be specified at runtime.

//    Performance tuning note:  For problems of substantial size,
//    preallocation of matrix memory is crucial for attaining good
//    performance. See the matrix chapter of the users manual for details.

//    We pass in nlocal as the "local" size of the matrix to force it
//    to have the same parallel layout as the vector created above.
//    */
//   if (!allocated_) {

//     ierr = MatCreate(PETSC_COMM_WORLD,&A_);CHKERRCONTINUE(ierr);
//     ierr = MatSetSizes(A_,nlocal,nlocal,n,n);CHKERRCONTINUE(ierr);
//     ierr = MatSetFromOptions(A_);CHKERRCONTINUE(ierr);
//     ierr = MatSetUp(A_);CHKERRCONTINUE(ierr);
//   } else {
//     // zero-out matrix
//     MatZeroEntries(A_);
//   }

//   /*
//    Assemble matrix.

//    The linear system is distributed across the processors by
//    chunks of contiguous rows, which correspond to contiguous
//    sections of the mesh on which the problem is discretized.
//    For matrix assembly, each processor contributes entries for
//    the part that it owns locally.
//    */

// #ifdef CPPUTILS_VERBOSE
//   out<<"  Zeroed-out sparse matrix entries... "<<timer<<endl;
// #endif

//   for (sparse_matrix_type::const_hash_iterator it = AA.map_.begin(); it !=
//   AA.map_.end(); ++it) {

//     // get subscripts
//     std::pair<size_t,size_t> subs = AA.unhash(it->first);
//     PetscInt row = subs.first;
//     PetscInt col = subs.second;
//     ierr = MatSetValues(A_, 1, &row, 1, &col, &it->second,
//     ADD_VALUES);CHKERRCONTINUE(ierr);
//   }

// #ifdef CPPUTILS_VERBOSE
//   out<<"  Filled sparse matrix... "<<timer<<endl;
// #endif

//   /* Assemble the matrix */
//   ierr = MatAssemblyBegin(A_,MAT_FINAL_ASSEMBLY);CHKERRCONTINUE(ierr);
//   ierr = MatAssemblyEnd(A_,MAT_FINAL_ASSEMBLY);CHKERRCONTINUE(ierr);

//   if (!allocated_) {
//     // set after the first MatAssemblyEnd
//     //        ierr = MatSetOption(A_, MAT_NEW_NONZERO_LOCATIONS,
//     PETSC_FALSE);CHKERRCONTINUE(ierr);
//     ierr = MatSetOption(A_, MAT_NEW_NONZERO_ALLOCATION_ERR,
//     PETSC_FALSE);CHKERRCONTINUE(ierr);
//   }

// #ifdef CPPUTILS_VERBOSE
//   out<<"  Sparse matrix assembled... "<<timer<<endl;
//   // view matrix
//   MatView(A_, PETSC_VIEWER_STDOUT_WORLD);

//   MatNorm(A_,NORM_FROBENIUS,&nrm);
//   out<<"  Norm of stiffness matrix... "<<nrm<<endl;
// #endif

//   /*
//    Set operators. Here the matrix that defines the linear system
//    also serves as the preconditioning matrix.
//    */
//   //    ierr =
//   KSPSetOperators(ksp,A_,A_,DIFFERENT_NONZERO_PATTERN);CHKERRCONTINUE(ierr);
//   ierr =
//   KSPSetOperators(ksp_,A_,A_,SAME_NONZERO_PATTERN);CHKERRCONTINUE(ierr);

//   //
//   //    /*
//   //     Set runtime options, e.g.,
//   //     -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
//   //     These options will override those specified above as long as
//   //     KSPSetFromOptions() is called _after_ any other customization
//   //     routines.
//   //     */
//   //    ierr = KSPSetFromOptions(ksp);CHKERRCONTINUE(ierr);

// #ifdef CPPUTILS_VERBOSE
//   out<<"  Solving system... "<<timer<<endl;
// #endif

//   /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    Solve the linear system
//    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
//   /*
//    Solve linear system
//    */
//   ierr = KSPSolve(ksp_,b,x_);CHKERRCONTINUE(ierr);

// #ifdef CPPUTILS_VERBOSE

//   /*
//    View solver info; we could instead use the option -ksp_view to
//    print this info to the screen at the conclusion of KSPSolve().
//    */
//   ierr = KSPView(ksp_,PETSC_VIEWER_STDOUT_WORLD);CHKERRCONTINUE(ierr);

//   int iter;
//   KSPGetIterationNumber(ksp_, &iter);
//   out<<"  System solved in "<<iter<<" iterations... "<<timer<<endl;
//   KSPConvergedReason reason;
//   ierr = KSPGetConvergedReason(ksp_,&reason);CHKERRCONTINUE(ierr);
//   if (reason < 0)
//     out<<"*** WARNING *** PETSc solver diverged with flag ";
//   else
//     out<<"*** INFO *** PETSc solver converged with flag ";

//   if (reason == KSP_CONVERGED_RTOL)
//     out<<"KSP_CONVERGED_RTOL"<<endl;
//   else if (reason == KSP_CONVERGED_ATOL)
//     out<<"KSP_CONVERGED_ATOL"<<endl;
//   else if (reason == KSP_CONVERGED_ITS)
//     out<<"KSP_CONVERGED_ITS"<<endl;
//   else if (reason == KSP_CONVERGED_CG_NEG_CURVE)
//     out<<"KSP_CONVERGED_CG_NEG_CURVE"<<endl;
//   else if (reason == KSP_CONVERGED_CG_CONSTRAINED)
//     out<<"KSP_CONVERGED_CG_CONSTRAINED"<<endl;
//   else if (reason == KSP_CONVERGED_STEP_LENGTH)
//     out<<"KSP_CONVERGED_STEP_LENGTH"<<endl;
//   else if (reason == KSP_CONVERGED_HAPPY_BREAKDOWN)
//     out<<"KSP_CONVERGED_HAPPY_BREAKDOWN"<<endl;
//   else if (reason == KSP_DIVERGED_NULL)
//     out<<"KSP_DIVERGED_NULL"<<endl;
//   else if (reason == KSP_DIVERGED_ITS)
//     out<<"KSP_DIVERGED_ITS"<<endl;
//   else if (reason == KSP_DIVERGED_DTOL)
//     out<<"KSP_DIVERGED_DTOL"<<endl;
//   else if (reason == KSP_DIVERGED_BREAKDOWN)
//     out<<"KSP_DIVERGED_BREAKDOWN"<<endl;
//   else if (reason == KSP_DIVERGED_BREAKDOWN_BICG)
//     out<<"KSP_DIVERGED_BREAKDOWN_BICG"<<endl;
//   else if (reason == KSP_DIVERGED_NONSYMMETRIC)
//     out<<"KSP_DIVERGED_NONSYMMETRIC"<<endl;
//   else if (reason == KSP_DIVERGED_INDEFINITE_PC)
//     out<<"KSP_DIVERGED_INDEFINITE_PC"<<endl;
//   else if (reason == KSP_DIVERGED_NAN)
//     out<<"KSP_DIVERGED_NAN"<<endl;
//   else if (reason == KSP_DIVERGED_INDEFINITE_MAT)
//     out<<"KSP_DIVERGED_INDEFINITE_MAT"<<endl;
//   else if (reason == KSP_CONVERGED_ITERATING)
//     out<<"KSP_CONVERGED_ITERATING"<<endl;

//   PetscReal rnorm;
//   ierr = KSPGetResidualNorm(ksp_,&rnorm);CHKERRCONTINUE(ierr);

//   out<<"PETSc last residual norm computed: "<<rnorm<<endl;

//   ierr = VecView(x_,PETSC_VIEWER_STDOUT_WORLD);CHKERRCONTINUE(ierr);

//   VecNorm(x_,NORM_2,&nrm);
//   out<<"  Norm of solution vector... "<<nrm<<endl;

// #endif

//   /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    Check solution and clean up
//    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

//   Vec x_all;

//   ierr = VecScatterCreateToAll(x_, &ctx, &x_all);CHKERRCONTINUE(ierr);
//   ierr =
//   VecScatterBegin(ctx,x_,x_all,INSERT_VALUES,SCATTER_FORWARD);CHKERRCONTINUE(ierr);
//   ierr =
//   VecScatterEnd(ctx,x_,x_all,INSERT_VALUES,SCATTER_FORWARD);CHKERRCONTINUE(ierr);

// #ifdef CPPUTILS_VERBOSE
//   out<<"  Solution vector scattered... "<<timer<<endl;
//   VecNorm(x_all,NORM_2,&nrm);
//   out<<"  Norm of scattered vector... "<<nrm<<endl;
//   //    ierr = VecView(x_all,PETSC_VIEWER_STDOUT_WORLD);CHKERRCONTINUE(ierr);
// #endif

//   /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    Get values from solution and store them in the object that will be
//    returned
//    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

//   sparse_vector_type xx(bb.size());

//   /* Set solution vector */
//   double zero = 0.;
//   double val;
//   for (sparse_vector_type::const_hash_iterator it = bb.map_.begin(); it !=
//   bb.map_.end(); ++it) {
//     int row = it->first;
//     ierr = VecGetValues(x_all, 1, &row, &val);
//     if (val != zero)
//       xx[row] = val;
//   }

// #ifdef CPPUTILS_VERBOSE
//   out<<"  Solution vector copied... "<<timer<<endl;
//   //    out<<"  Norm of copied solution vector... "<<norm(xx,
//   Insert_t)<<endl;
// #endif

//   /*
//    Free work space.  All PETSc objects should be destroyed when they
//    are no longer needed.
//    */
//   ierr = VecDestroy(&b);CHKERRCONTINUE(ierr);
//   //    ierr = KSPDestroy(&ksp);CHKERRCONTINUE(ierr);

//   // set allocated flag
//   if (!allocated_) {
//     allocated_ = true;
//   }

// #ifdef CPPUTILS_VERBOSE
//   out<<"  Temporary data structures destroyed... "<<timer<<endl;
// #endif

//   return xx;
// }

// SolverPETSc::sparse_vector_type SolverPETSc::operator()(const
// SolverPETSc::sparse_matrix_type& KKpf, const SolverPETSc::sparse_matrix_type&
// KKpp, const SolverPETSc::sparse_vector_type& UUp) {

// #ifdef CPPUTILS_VERBOSE
//   // parallel output stream
//   Output_stream out;
//   // timer
//   cpputils::ctimer timer;
//   out<<"Inside SolverPETSc::operator()(const sparse_matrix&, const
//   sparse_matrix&, const sparse_vector&). "<<timer<<endl;
// #endif

//   Mat Kpf, Kpp;
//   Vec Up, Pf, Pp;

//   PetscErrorCode ierr =
//   MatCreate(PETSC_COMM_WORLD,&Kpf);CHKERRCONTINUE(ierr);
//   ierr = MatCreate(PETSC_COMM_WORLD,&Kpp);CHKERRCONTINUE(ierr);

// #ifdef CPPUTILS_VERBOSE
//   out<<"  Allocating memory... "<<timer<<endl;
// #endif

//   ierr = MatSetFromOptions(Kpf);CHKERRCONTINUE(ierr);
//   ierr = MatSetFromOptions(Kpp);CHKERRCONTINUE(ierr);

//   ierr = MatSetSizes(Kpf,PETSC_DECIDE,PETSC_DECIDE, KKpf.rows(),
//   KKpf.columns());CHKERRCONTINUE(ierr);
//   ierr = MatSetSizes(Kpp,PETSC_DECIDE,PETSC_DECIDE, KKpp.rows(),
//   KKpp.columns());CHKERRCONTINUE(ierr);

//   // get number of non-zeros in both diagonal and non-diagonal portions of
//   the matrix

//   std::pair<size_t,size_t> Kpf_nz = KKpf.non_zero_off_diagonal();
//   std::pair<size_t,size_t> Kpp_nz = KKpp.non_zero_off_diagonal();

//   ierr = MatMPIAIJSetPreallocation(Kpf, Kpf_nz.first, PETSC_NULL,
//   Kpf_nz.second, PETSC_NULL); CHKERRCONTINUE(ierr);
//   ierr = MatMPIAIJSetPreallocation(Kpp, Kpp_nz.first, PETSC_NULL,
//   Kpp_nz.second, PETSC_NULL); CHKERRCONTINUE(ierr);
//   ierr = MatSeqAIJSetPreallocation(Kpf, Kpf_nz.first, PETSC_NULL);
//   CHKERRCONTINUE(ierr);
//   ierr = MatSeqAIJSetPreallocation(Kpp, Kpp_nz.first, PETSC_NULL);
//   CHKERRCONTINUE(ierr);

//   for (sparse_matrix_type::const_hash_iterator it = KKpf.map_.begin(); it !=
//   KKpf.map_.end(); ++it) {

//     // get subscripts
//     std::pair<size_t,size_t> subs = KKpf.unhash(it->first);
//     int row = subs.first;
//     int col = subs.second;
//     ierr = MatSetValues(Kpf, 1, &row, 1, &col, &it->second,
//     ADD_VALUES);CHKERRCONTINUE(ierr);
//   }

//   for (sparse_matrix_type::const_hash_iterator it = KKpp.map_.begin(); it !=
//   KKpp.map_.end(); ++it) {

//     // get subscripts
//     std::pair<size_t,size_t> subs = KKpp.unhash(it->first);
//     int row = subs.first;
//     int col = subs.second;
//     ierr = MatSetValues(Kpf, 1, &row, 1, &col, &it->second,
//     ADD_VALUES);CHKERRCONTINUE(ierr);
//   }

// #ifdef CPPUTILS_VERBOSE
//   out<<"  Filled sparse matrices... "<<timer<<endl;
// #endif

//   /*
//    Assemble matrix, using the 2-step process:
//    MatAssemblyBegin(), MatAssemblyEnd()
//    Computations can be done while messages are in transition
//    by placing code between these two statements.
//    */
//   ierr = MatAssemblyBegin(Kpf,MAT_FINAL_ASSEMBLY);CHKERRCONTINUE(ierr);
//   ierr = MatAssemblyBegin(Kpp,MAT_FINAL_ASSEMBLY);CHKERRCONTINUE(ierr);
//   ierr = MatAssemblyEnd(Kpf,MAT_FINAL_ASSEMBLY);CHKERRCONTINUE(ierr);
//   ierr = MatAssemblyEnd(Kpp,MAT_FINAL_ASSEMBLY);CHKERRCONTINUE(ierr);

// #ifdef CPPUTILS_VERBOSE
//   out<<"  Sparse matrices assembled... "<<timer<<endl;
// #endif

//   ierr = VecCreate(PETSC_COMM_WORLD,&Up);CHKERRCONTINUE(ierr);
//   ierr = VecSetSizes(Up,PETSC_DECIDE, UUp.size());CHKERRCONTINUE(ierr);
//   ierr = VecSetFromOptions(Up);CHKERRCONTINUE(ierr);

//   ierr = VecCreate(PETSC_COMM_WORLD,&Pf);CHKERRCONTINUE(ierr);
//   ierr = VecSetSizes(Pf,PETSC_DECIDE, KKpf.rows());CHKERRCONTINUE(ierr);
//   ierr = VecSetFromOptions(Pf);CHKERRCONTINUE(ierr);
//   ierr = VecDuplicate(Pf,&Pp);CHKERRCONTINUE(ierr);

// #ifdef CPPUTILS_VERBOSE
//   out<<"  Vectors created... "<<timer<<endl;
// #endif

//   /* Set hight-hand-side vector */
//   for (sparse_vector_type::const_hash_iterator it = UUp.map_.begin(); it !=
//   UUp.map_.end(); ++it) {
//     int row = it->first;
//     ierr = VecSetValues(Up, 1, &row, &it->second, ADD_VALUES);
//   }

//   /*
//    Assemble vector, using the 2-step process:
//    VecAssemblyBegin(), VecAssemblyEnd()
//    Computations can be done while messages are in transition
//    by placing code between these two statements.
//    */
//   ierr = VecAssemblyBegin(Up);CHKERRCONTINUE(ierr);
//   ierr = VecAssemblyEnd(Up);CHKERRCONTINUE(ierr);

//   // add Kpf*Uf
//   ierr = MatMult(Kpf, x_, Pf);

//   // add Kpp*Up
//   ierr = MatMultAdd(Kpp, Up, Pf, Pp);

// #ifdef CPPUTILS_VERBOSE
//   out<<"  Matrices multiplied... "<<timer<<endl;
// #endif

//   VecScatter ctx;
//   Vec Pp_all;

//   ierr = VecScatterCreateToAll(Pp, &ctx, &Pp_all);CHKERRCONTINUE(ierr);
//   ierr =
//   VecScatterBegin(ctx,Pp,Pp_all,INSERT_VALUES,SCATTER_FORWARD);CHKERRCONTINUE(ierr);
//   ierr =
//   VecScatterEnd(ctx,Pp,Pp_all,INSERT_VALUES,SCATTER_FORWARD);CHKERRCONTINUE(ierr);

// #ifdef CPPUTILS_VERBOSE
//   out<<"  Vector scattered... "<<timer<<endl;
// #endif

//   /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    Get values from solution and store them in the object that will be
//    returned
//    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

//   sparse_vector_type pp(KKpf.rows());

//   // get reaction vector
//   for (int i=0; i<KKpf.rows(); ++i) {

//     PetscScalar v;
//     ierr = VecGetValues(Pp_all, 1, &i, &v);
//     if (v != PetscScalar())
//       pp[i] = v;
//   }

// #ifdef CPPUTILS_VERBOSE
//   out<<"  Vector copied... "<<timer<<endl;
// #endif

//   ierr = MatDestroy(&Kpf);CHKERRCONTINUE(ierr);
//   ierr = MatDestroy(&Kpp);CHKERRCONTINUE(ierr);
//   ierr = VecDestroy(&Up);CHKERRCONTINUE(ierr);
//   ierr = VecDestroy(&Pf);CHKERRCONTINUE(ierr);
//   ierr = VecDestroy(&Pp);CHKERRCONTINUE(ierr);
//   ierr = VecDestroy(&Pp_all);CHKERRCONTINUE(ierr);

// #ifdef CPPUTILS_VERBOSE
//   out<<"  Temporary data structures destroyed... "<<timer<<endl;
// #endif

//   return pp;
// }

// SolverPETSc::sparse_vector_type SolverPETSc::operator()(const
// SolverPETSc::sparse_vector_type& aa, const SolverPETSc::sparse_vector_type&
// bb) {

//   assert(aa.size() == bb.size());

// #ifdef CPPUTILS_VERBOSE
//   // parallel output stream
//   Output_stream out;
//   // timer
//   cpputils::ctimer timer;
//   out<<"Inside SolverPETSc::operator()(const sparse_vector&, const
//   sparse_vector&). "<<timer<<endl;
// #endif

//   Vec r;

//   PetscErrorCode ierr = VecCreate(PETSC_COMM_WORLD,&r);CHKERRCONTINUE(ierr);
//   ierr = VecSetSizes(r,PETSC_DECIDE, aa.size());CHKERRCONTINUE(ierr);
//   ierr = VecSetFromOptions(r);CHKERRCONTINUE(ierr);

// #ifdef CPPUTILS_VERBOSE
//   out<<"  Vectors created... "<<timer<<endl;
// #endif

//   // set values
//   for (sparse_vector_type::const_hash_iterator it = aa.map_.begin(); it !=
//   aa.map_.end(); ++it) {
//     int row = it->first;
//     ierr = VecSetValues(r, 1, &row, &it->second, ADD_VALUES);
//   }
//   for (sparse_vector_type::const_hash_iterator it = bb.map_.begin(); it !=
//   bb.map_.end(); ++it) {
//     int row = it->first;
//     ierr = VecSetValues(r, 1, &row, &it->second, ADD_VALUES);
//   }

//   /*
//    Assemble vector, using the 2-step process:
//    VecAssemblyBegin(), VecAssemblyEnd()
//    Computations can be done while messages are in transition
//    by placing code between these two statements.
//    */
//   ierr = VecAssemblyBegin(r);CHKERRCONTINUE(ierr);
//   ierr = VecAssemblyEnd(r);CHKERRCONTINUE(ierr);

//   sparse_vector_type rr(aa.size());

//   for (sparse_vector_type::const_hash_iterator it = aa.map_.begin(); it !=
//   aa.map_.end(); ++it) {
//     int row = it->first;
//     ierr = VecGetValues(r, 1, &row, &rr[row]);
//   }
//   for (sparse_vector_type::const_hash_iterator it = bb.map_.begin(); it !=
//   bb.map_.end(); ++it) {
//     int row = it->first;
//     ierr = VecGetValues(r, 1, &row, &rr[row]);
//   }

// #ifdef CPPUTILS_VERBOSE
//   out<<"  Vector copied... "<<timer<<endl;
// #endif

//   ierr = VecDestroy(&r);CHKERRCONTINUE(ierr);

// #ifdef CPPUTILS_VERBOSE
//   out<<"  Temporary data structures destroyed... "<<timer<<endl;
// #endif

//   return rr;
// }

// SolverPETSc::value_type SolverPETSc::norm(const
// SolverPETSc::sparse_matrix_type& aa, Element_insertion_type flag) {

// #ifdef CPPUTILS_VERBOSE
//   // parallel output stream
//   Output_stream out;
//   // timer
//   cpputils::ctimer timer;
//   out<<"Inside SolverPETSc::norm(const sparse_matrix&). "<<timer<<endl;
// #endif

//   Mat r;

//   PetscErrorCode ierr = MatCreate(PETSC_COMM_WORLD,&r);CHKERRCONTINUE(ierr);
//   ierr = MatSetSizes(r,PETSC_DECIDE,PETSC_DECIDE, aa.rows(),
//   aa.columns());CHKERRCONTINUE(ierr);
//   ierr = MatSetFromOptions(r);CHKERRCONTINUE(ierr);

// #ifdef CPPUTILS_VERBOSE
//   out<<"  Matrix created... "<<timer<<endl;
// #endif

//   // set values
//   for (sparse_vector_type::const_hash_iterator it = aa.map_.begin(); it !=
//   aa.map_.end(); ++it) {
//     // get subscripts
//     std::pair<size_t,size_t> subs = aa.unhash(it->first);
//     int row = subs.first;
//     int col = subs.second;

//     if (flag == Add_t)
//       ierr = MatSetValues(r, 1, &row, 1, &col, &it->second, ADD_VALUES);
//     else if (flag == Insert_t)
//       ierr = MatSetValues(r, 1, &row, 1, &col, &it->second, INSERT_VALUES);
//     CHKERRCONTINUE(ierr);
//   }

// #ifdef CPPUTILS_VERBOSE
//   out<<"  Matrix filled..."<<timer<<endl;
// #endif

//   /*
//    Assemble vector, using the 2-step process:
//    VecAssemblyBegin(), VecAssemblyEnd()
//    Computations can be done while messages are in transition
//    by placing code between these two statements.
//    */
//   ierr = MatAssemblyBegin(r,MAT_FINAL_ASSEMBLY);CHKERRCONTINUE(ierr);
//   ierr = MatAssemblyEnd(r,MAT_FINAL_ASSEMBLY);CHKERRCONTINUE(ierr);

//   value_type nrm;

//   MatNorm(r,NORM_FROBENIUS,&nrm);

// #ifdef CPPUTILS_VERBOSE
//   out<<"  Norm computed... "<<timer<<endl;
// #endif

//   ierr = MatDestroy(&r);CHKERRCONTINUE(ierr);

// #ifdef CPPUTILS_VERBOSE
//   out<<"  Temporary data structures destroyed... "<<timer<<endl;
// #endif

//   return nrm;
// }

// SolverPETSc::value_type SolverPETSc::norm(const
// SolverPETSc::sparse_vector_type& aa, Element_insertion_type flag) {

// #ifdef CPPUTILS_VERBOSE
//   // parallel output stream
//   Output_stream out;
//   // timer
//   cpputils::ctimer timer;
//   out<<"Inside SolverPETSc::norm(const sparse_vector&). "<<timer<<endl;
// #endif

//   Vec r;

//   PetscErrorCode ierr = VecCreate(PETSC_COMM_WORLD,&r);CHKERRCONTINUE(ierr);
//   ierr = VecSetSizes(r,PETSC_DECIDE, aa.size());CHKERRCONTINUE(ierr);
//   ierr = VecSetFromOptions(r);CHKERRCONTINUE(ierr);

// #ifdef CPPUTILS_VERBOSE
//   out<<"  Vector created... "<<timer<<endl;
// #endif

//   // set values
//   for (sparse_vector_type::const_hash_iterator it = aa.map_.begin(); it !=
//   aa.map_.end(); ++it) {
//     int row = it->first;
//     if (flag == Add_t)
//       ierr = VecSetValues(r, 1, &row, &it->second, ADD_VALUES);
//     else if (flag == Insert_t)
//       ierr = VecSetValues(r, 1, &row, &it->second, INSERT_VALUES);
//     CHKERRCONTINUE(ierr);
//   }

// #ifdef CPPUTILS_VERBOSE
//   out<<"  Vector filled..."<<timer<<endl;
// #endif

//   /*
//    Assemble vector, using the 2-step process:
//    VecAssemblyBegin(), VecAssemblyEnd()
//    Computations can be done while messages are in transition
//    by placing code between these two statements.
//    */
//   ierr = VecAssemblyBegin(r);CHKERRCONTINUE(ierr);
//   ierr = VecAssemblyEnd(r);CHKERRCONTINUE(ierr);

//   value_type nrm;

//   VecNorm(r,NORM_2,&nrm);

// #ifdef CPPUTILS_VERBOSE
//   out<<"  Norm computed... "<<timer<<endl;
// #endif

//   ierr = VecDestroy(&r);CHKERRCONTINUE(ierr);

// #ifdef CPPUTILS_VERBOSE
//   out<<"  Temporary data structures destroyed... "<<timer<<endl;
// #endif

//   return nrm;

// }

// //
// ///*
// -------------------------------------------------------------------------- */
// //SolverMumps::SolverMumps(SparseMatrix & matrix,
// //                         const ID & id,
// //                         const MemoryID & memory_id) :
// //Solver(matrix, id, memory_id), is_mumps_data_initialized(false),
// rhs_is_local(true) {
// //  AKANTU_DEBUG_IN();
// //
// //#ifdef AKANTU_USE_MPI
// //  parallel_method = SolverMumpsOptions::_fully_distributed;
// //#else //AKANTU_USE_MPI
// //  parallel_method = SolverMumpsOptions::_master_slave_distributed;
// //#endif //AKANTU_USE_MPI
// //
// //  CommunicatorEventHandler & comm_event_handler = *this;
// //
// //  communicator.registerEventHandler(comm_event_handler);
// //
// //  AKANTU_DEBUG_OUT();
// //}
// //
// ///*
// -------------------------------------------------------------------------- */
// //SolverMumps::~SolverMumps() {
// //  AKANTU_DEBUG_IN();
// //
// //  AKANTU_DEBUG_OUT();
// //}
// //
// ///*
// -------------------------------------------------------------------------- */
// //void SolverMumps::destroyMumpsData() {
// //  AKANTU_DEBUG_IN();
// //
// //  if(is_mumps_data_initialized) {
// //    mumps_data.job = _smj_destroy; // destroy
// //    dmumps_c(&mumps_data);
// //    is_mumps_data_initialized = false;
// //  }
// //
// //  AKANTU_DEBUG_OUT();
// //}
// //
// ///*
// -------------------------------------------------------------------------- */
// //void SolverMumps::onCommunicatorFinalize(const StaticCommunicator & comm) {
// //  AKANTU_DEBUG_IN();
// //
// //  try{
// //    const StaticCommunicatorMPI & comm_mpi =
// //    dynamic_cast<const StaticCommunicatorMPI
// &>(comm.getRealStaticCommunicator());
// //    if(mumps_data.comm_fortran ==
// MPI_Comm_c2f(comm_mpi.getMPICommunicator()))
// //      destroyMumpsData();
// //  } catch(...) {}
// //
// //  AKANTU_DEBUG_OUT();
// //}
// //
// ///*
// -------------------------------------------------------------------------- */
// //void SolverMumps::initMumpsData(SolverMumpsOptions::ParallelMethod
// parallel_method) {
// //  switch(parallel_method) {
// //    case SolverMumpsOptions::_fully_distributed:
// //      icntl(18) = 3; //fully distributed
// //      icntl(28) = 0; //automatic choice
// //
// //      mumps_data.nz_loc  = matrix->getNbNonZero();
// //      mumps_data.irn_loc = matrix->getIRN().values;
// //      mumps_data.jcn_loc = matrix->getJCN().values;
// //      break;
// //    case SolverMumpsOptions::_master_slave_distributed:
// //      if(prank == 0) {
// //        mumps_data.nz  = matrix->getNbNonZero();
// //        mumps_data.irn = matrix->getIRN().values;
// //        mumps_data.jcn = matrix->getJCN().values;
// //      } else {
// //        mumps_data.nz  = 0;
// //        mumps_data.irn = NULL;
// //        mumps_data.jcn = NULL;
// //
// //        icntl(18) = 0; //centralized
// //        icntl(28) = 0; //sequential analysis
// //      }
// //      break;
// //  }
// //}
// //
// ///*
// -------------------------------------------------------------------------- */
// //void SolverMumps::initialize(SolverOptions & options) {
// //  AKANTU_DEBUG_IN();
// //
// //  mumps_data.par = 1;
// //
// //  if(SolverMumpsOptions * opt = dynamic_cast<SolverMumpsOptions
// *>(&options)) {
// //    if(opt->parallel_method ==
// SolverMumpsOptions::_master_slave_distributed) {
// //      mumps_data.par = 0;
// //    }
// //  }
// //
// //  mumps_data.sym = 2 * (matrix->getSparseMatrixType() == _symmetric);
// //  prank = communicator.whoAmI();
// //#ifdef AKANTU_USE_MPI
// //  mumps_data.comm_fortran = MPI_Comm_c2f(dynamic_cast<const
// StaticCommunicatorMPI
// &>(communicator.getRealStaticCommunicator()).getMPICommunicator());
// //#endif
// //
// //  if(AKANTU_DEBUG_TEST(dblTrace)) {
// //    icntl(1) = 2;
// //    icntl(2) = 2;
// //    icntl(3) = 2;
// //    icntl(4) = 4;
// //  }
// //
// //  mumps_data.job = _smj_initialize; //initialize
// //  dmumps_c(&mumps_data);
// //  is_mumps_data_initialized = true;
// //
// //  /*
// ------------------------------------------------------------------------ */
// //  UInt size = matrix->size();
// //
// //  if(prank == 0) {
// //    std::stringstream sstr_rhs; sstr_rhs << id << ":rhs";
// //    rhs = &(alloc<Real>(sstr_rhs.str(), size, 1, REAL_INIT_VALUE));
// //  } else {
// //    rhs = NULL;
// //  }
// //
// //  /// No outputs
// //  icntl(1) = 0;
// //  icntl(2) = 0;
// //  icntl(3) = 0;
// //  icntl(4) = 0;
// //  mumps_data.nz_alloc = 0;
// //
// //  if (AKANTU_DEBUG_TEST(dblDump)) icntl(4) = 4;
// //
// //  mumps_data.n   = size;
// //
// //  if(AKANTU_DEBUG_TEST(dblDump)) {
// //    strcpy(mumps_data.write_problem, "mumps_matrix.mtx");
// //  }
// //
// //  /*
// ------------------------------------------------------------------------ */
// //  // Default Scaling
// //  icntl(8) = 77;
// //
// //  icntl(5) = 0; // Assembled matrix
// //
// //  SolverMumpsOptions * opt = dynamic_cast<SolverMumpsOptions *>(&options);
// //  if(opt)
// //    parallel_method = opt->parallel_method;
// //
// //  initMumpsData(parallel_method);
// //
// //  mumps_data.job = _smj_analyze; //analyze
// //  dmumps_c(&mumps_data);
// //
// //  AKANTU_DEBUG_OUT();
// //}
// //
// ///*
// -------------------------------------------------------------------------- */
// //void SolverMumps::setRHS(Array<Real> & rhs) {
// //  if(prank == 0) {
// //    matrix->getDOFSynchronizer().gather(rhs, 0, this->rhs);
// //  } else {
// //    matrix->getDOFSynchronizer().gather(rhs, 0);
// //  }
// //}
// //
// ///*
// -------------------------------------------------------------------------- */
// //void SolverMumps::solve() {
// //  AKANTU_DEBUG_IN();
// //
// //  if(parallel_method == SolverMumpsOptions::_fully_distributed)
// //    mumps_data.a_loc  = matrix->getA().values;
// //  else
// //    if(prank == 0) {
// //      mumps_data.a  = matrix->getA().values;
// //    }
// //
// //  if(prank == 0) {
// //    mumps_data.rhs = rhs->values;
// //  }
// //
// //  /// Default centralized dense second member
// //  icntl(20) = 0;
// //  icntl(21) = 0;
// //
// //  mumps_data.job = _smj_factorize_solve; //solve
// //  dmumps_c(&mumps_data);
// //
// //  AKANTU_DEBUG_ASSERT(info(1) != -10, "Singular matrix");
// //  AKANTU_DEBUG_ASSERT(info(1) == 0,
// //                      "Error in mumps during solve process, check mumps
// user guide INFO(1) ="
// //                      << info(1));
// //
// //  AKANTU_DEBUG_OUT();
// //}
// //
// ///*
// -------------------------------------------------------------------------- */
// //void SolverMumps::solve(Array<Real> & solution) {
// //  AKANTU_DEBUG_IN();
// //
// //  solve();
// //
// //  if(prank == 0) {
// //    matrix->getDOFSynchronizer().scatter(solution, 0, this->rhs);
// //  } else {
// //    matrix->getDOFSynchronizer().scatter(solution, 0);
// //  }
// //
// //  AKANTU_DEBUG_OUT();
// //}

} // namespace akantu
