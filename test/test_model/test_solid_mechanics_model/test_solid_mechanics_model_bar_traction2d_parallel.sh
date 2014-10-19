#!/bin/bash


rm energy_bar_2d_para.csv

mpirun -np 2 ./test_solid_mechanics_model_bar_traction2d_parallel
ret=$?
if [ $ret -eq 0 ]
then
    ./test_cst_energy.pl energy_bar_2d_para.csv 1e-2
else
    return $ret
fi
