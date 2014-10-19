#!/bin/bash

rm energy_2d_pbc.csv
./test_solid_mechanics_model_bar_traction2d_structured_pbc
ret=$?
if [ $ret -eq 0 ]
then
    ./test_cst_energy.pl energy_2d_pbc.csv 1e-4
else
    return $ret
fi