#!/bin/bash

./test_igfem_element_orientation
mpirun -np 10 ./test_igfem_element_orientation
