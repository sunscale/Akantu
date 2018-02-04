#!/bin/bash

set -o errexit
set -o pipefail

show_help() {
    cat << EOF
Usage: ${0##*/} -n NAME -e EXECUTABLE [-p MPI_WRAPPER] [-s SCRIPT_FILE]
          [-r REFERENCE_FILE] [-w WORKING_DIR] [ARGS]
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
    -E ENVIRONMENT_FILE File to source before running tests
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
    echo "Run $*"
    (($* | tee "${name}${sout}") 3>&1 1>&2 2>&3 | tee "${name}${serr}") 3>&1 1>&2 2>&3

    lastout="${name}${sout}"
}

name=
executable=
parallel=
postprocess_script=
reference=
working_dir=
envi=
parallel_processes="2"

while getopts ":n:e:E:p:N:s:r:w:h" opt; do
    case "$opt" in
        n)  name="$OPTARG"
            ;;
        e)  executable="$OPTARG"
            ;;
        p)  parallel="$OPTARG"
            ;;
        N)  parallel_processes="$OPTARG"
            ;;
        s)  postprocess_script="$OPTARG"
            ;;
        r)  reference="$OPTARG"
            ;;
        w)  working_dir="$OPTARG"
            ;;
        E)  envi="$OPTARG"
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

shift $(( $OPTIND - 1 ))
_args=$*

if [ -n "${envi}" ]; then
    source ${envi}
fi

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
    full_redirect 0 ${name} "${executable} ${_args}"
else
  #for i in ${parallel_processes}; do
  i=${parallel_processes}
  echo "Executing the test ${name} for ${i} procs"
  full_redirect $i ${name}_$i "${parallel}  ${i} ${executable} ${_args}"
  #done
fi

if [ -n "${postprocess_script}" ]; then
  echo "Executing the test ${name} post-processing"
  full_redirect 0 ${name}_pp ./${postprocess_script}
fi

if [ -n "${reference}" ]; then
   echo "Comparing last generated output to the reference file"
   diff ${lastout} ${reference}
fi
