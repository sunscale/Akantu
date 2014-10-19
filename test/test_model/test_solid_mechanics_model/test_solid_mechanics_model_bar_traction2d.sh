#!/bin/bash

rm energy_bar_2d.csv
./test_solid_mechanics_model_bar_traction2d
ret=$?
if [ $ret -eq 0 ]
then
    ./test_cst_energy.pl energy_bar_2d.csv 1e-2
else
    return $ret
fi