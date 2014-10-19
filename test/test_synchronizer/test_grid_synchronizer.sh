#!/bin/bash

./test_grid_synchronizer
cp neighbors_ref_0 neighbors_ref

mpirun -np 8 ./test_grid_synchronizer
mpirun -np 8 ./test_grid_synchronizer_check_neighbors
