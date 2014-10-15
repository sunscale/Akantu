#! /bin/bash

rm -fr paraview
mkdir paraview

rm -fr displacement
mkdir displacement

./test_cohesive_parallel_extrinsic_tetrahedron_displacement
mpirun -np 8 ./test_cohesive_parallel_extrinsic_tetrahedron_displacement
