#!/bin/bash

./$1 > $1.lastout 2> $1.errout

ret=$?

if [ $ret -eq 0 ]
then
    diff $1.lastout $2 && echo "Test passed!!!"
else
    echo "Test Failed!!"
    exit $ret
fi
