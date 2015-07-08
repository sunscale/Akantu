#!/bin/bash

./$1 > $1.lastout 2> $1.errout

ret=$?

if [ $ret -eq 0 ]
then
    cat $1.lastout | tr '\n' '\r' | sed "s/--------------------------------------------------------------------------\r\[\[[0-9]\+,[0-9]\+\],[0-9]\+\]: A high-performance Open MPI point-to-point messaging module\rwas unable to find any relevant network interfaces:\r\rModule: OpenFabrics (openib)\r  Host: scratch\r\rAnother transport will be used instead, although this may result in\rlower performance.\r--------------------------------------------------------------------------\r//" | tr '\r' '\n' > $1.lastout.replacement
    cp $1.lastout.replacement $1.lastout
    diff -q $1.lastout $2 && echo "Test passed!!!"
else
    echo "Test Failed!!"
    exit $ret
fi
