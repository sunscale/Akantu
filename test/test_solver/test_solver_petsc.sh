#! /bin/bash

# choose solver Mumps through the PETSc interface
./test_solver_petsc -ksp_type preonly -pc_type lu -pc_factor_mat_solver_package mumps -mat_mumps_icntl_7 2


