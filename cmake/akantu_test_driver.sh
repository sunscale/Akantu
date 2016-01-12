#!/bin/bash

set -o errexit

show_help() {
    cat << EOF
Usage: ${0##*/} -n NAME -e EXECUTABLE [-p MPI_WRAPPER] [-s SCRIPT_FILE]
          [-r REFERENCE_FILE] [-w WORKING_DIR]
Execute the test in the good configuration according to the options given

    -e EXECUTABLE     Main executable of the test
    -n NAME           Name of the test
    -p MPI_WRAPPER    Executes the test for multiple parallel configuration
    -s SCRIPT_FILE    Script to execute after the execution of the test to
                      postprocess the results
    -r REFERENCE_FILE Reference file to compare with if the name of the file
                      contains a <nb_proc> this will be used for the different
                      configuration when -p is given
    -w WORKING_DIR    The directory in which to execute the test
    -h                Print this helps
EOF
}

full_redirect() {
    local nproc=$1
    shift
    local name=$1
    shift

    local sout=".lastout"
    local serr=".lasterr"
    if [ ${nproc} -ne 0 ]; then
       sout="-${nproc}${sout}"
       serr="-${nproc}${serr}"
    fi
    (($* | tee "${name}${sout}") 3>&1 1>&2 2>&3 | tee "${name}${serr}") 3>&1 1>&2 2>&3

    res=$?
    if [ ! $res -eq 0 ]; then
        exit $res
    fi

    lastout="${name}${sout}"
}

name=
executable=
parallel=
postprocess_script=
reference=
working_dir=

while getopts ":n:e:p:s:r:w:h" opt; do
    case "$opt" in
        n)  name="$OPTARG"
            ;;
        e)  executable="$OPTARG"
            ;;
        p)  parallel="$OPTARG"
            ;;
        s)  postprocess_script="$OPTARG"
            ;;
        r)  reference="$OPTARG"
            ;;
        w)  working_dir="$OPTARG"
            ;;
        h)
            show_help
            exit 0
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            show_help
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            show_help
            exit 1
            ;;
    esac
done

if [ -z "${name}" -o -z "${executable}" ]; then
    echo "Missing executable or name"
    show_help
    exit 1
fi

if [ -n "${working_dir}" ]; then
    current_directory=$PWD
    echo "Entering directory ${working_dir}"
    cd "${working_dir}"
fi

if [ -z "${parallel}" ]; then
    echo "Executing the test ${name}"
    full_redirect 0 ${name} "./${executable}"
else
    for i in ${parallel_processes}; do
        echo "Executing the test ${name} for ${i} procs"
        full_redirect $i ${name}_$i "${parallel_processes}  ${i} ./${executable}"
    done
fi

if [ -n "${postprocess_script}" ]; then
    echo "Executing the test ${name} post-processing"
    full_redirect 0 ${name}_pp ./${postprocess_script}
fi

if [ -n "${reference}" ]; then
   echo "Comparing last generated output to the reference file"
   diff ${lastout} ${reference}
fi

