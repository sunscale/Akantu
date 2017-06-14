#include <stdio.h>

#if !defined(MUMPS_SEQ)
#  include <mpi.h>
#endif

#define JOB_INIT -1
#define JOB_END -2
#define JOB_COMPLETE 6
#define USE_COMM_WORLD -987654

#define icntl(n) id.icntl[n - 1]

int main(int argc, char **argv) {
  int n = 2;
  int nz = 2;

  int irn[2] = {1, 2};
  int jcn[2] = {1, 2};
  Real a[2];
  Real rhs[2];

#if !defined(MUMPS_SEQ)
  MPI_Init(&argc, &argv);
#endif

  rhs[0] = 1.0; rhs[1]=4.0;
  a[0] = 1.0; a[1] = 2.0;

  id.job = JOB_INIT;
  id.par = 1;
  id.sym = 0;

#if !defined(MUMPS_SEQ)
  id.comm_fortran = USE_COMM_WORLD;
#endif

  mumps_c(&id);

  // Default Scaling
  icntl(8) = 77;

  // Assembled matrix
  icntl(5) = 0;

  /// Default centralized dense second member
  icntl(20) = 0;
  icntl(21) = 0;

  // automatic choice for analysis analysis
  icntl(28) = 0;

  // fully distributed
  icntl(18) = 3;

  id.n = n;

  id.nz_loc = nz;
  id.irn_loc = irn;
  id.jcn_loc = jcn;

  id.a_loc = a;
  id.rhs = rhs;

  icntl(1) = -1;
  icntl(2) = -1;
  icntl(3) = -1;
  icntl(4) = 0;


  id.job = JOB_COMPLETE;
  mumps_c(&id);

  id.job=JOB_END;
  mumps_c(&id);

  printf("Solution is : (%8.2f %8.2f)\n", rhs[0], rhs[1]);

  return 0;
}
