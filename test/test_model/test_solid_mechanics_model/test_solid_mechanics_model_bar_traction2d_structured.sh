#!/bin/bash

rm energy.csv

./test_solid_mechanics_model_bar_traction2d_structured
ret=$?
if [ $ret -eq 0 ]
then
    ./test_cst_energy.pl energy_bar_2d_structured.csv 1e-3
else
    return $ret
fi
