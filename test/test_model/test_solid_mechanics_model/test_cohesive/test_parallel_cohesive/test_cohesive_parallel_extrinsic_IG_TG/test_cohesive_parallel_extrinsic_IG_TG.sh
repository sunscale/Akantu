#!/bin/bash

rm -r paraview
mkdir paraview

mpirun -np 8 ./test_cohesive_parallel_extrinsic_IG_TG
