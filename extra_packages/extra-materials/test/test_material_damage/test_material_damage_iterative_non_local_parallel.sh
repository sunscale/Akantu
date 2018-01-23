#! /bin/bash

rm -fr paraview
mkdir paraview

rm -fr displacement
mkdir displacement

./test_material_damage_iterative_non_local_parallel
mpirun -np 4 ./test_material_damage_iterative_non_local_parallel
