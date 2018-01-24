#!/bin/bash

status=0

PARALLEL_LEVEL=( 2 4 )
pwd
for (( p=0; p<$PARALLEL_LEVEL; ++p)); do
    PRANK=${PARALLEL_LEVEL[${p}]}
    
    for (( i=0; i<$PRANK; ++i)); do
	echo "Checking output_from_proc_"${i}"_out_of_"${PRANK}".out"
	FILE_NAME="output_from_proc_"${i}"_out_of_"${PRANK}".out"
	diff "output_dir/"${FILE_NAME} "output_dir_verified/"${FILE_NAME}
	status=$(($? + ${status}))
    done
done

exit ${status}
