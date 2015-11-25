#!/bin/bash

./test_pair_computation
mpirun -np 4 ./test_pair_computation
