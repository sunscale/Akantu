#!/bin/bash

./test_assembling_K_cohe_elements -ksp_type preonly -pc_type lu -pc_factor_mat_solver_package mumps -mat_mumps_icntl_7 2 -ksp_initial_guess_nonzero 0
