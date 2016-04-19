#! /bin/bash

mpirun -np 4 ./test_solver_petsc -ksp_type preonly -pc_type lu -pc_factor_mat_solver_package mumps -mat_mumps_icntl_7 2
