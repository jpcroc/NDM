#!/bin/bash -f
# **************************************************************************************************************************
#   Set environment
# **************************************************************************************************************************
# --------------------------------------- necessary environment to run or build ndm ----------------------------------------

# On Gatsby installed lammps 2021
module load lammps/29Sep2021-u2

# On Gatsby locally installed LAMMPS
#module load compiler
#module load mpi
#module load mkl

# ---------------------------------------- more environment, only relevant for build ---------------------------------------
export NDM_ROODIR=/volatile/catB/jd270899/git_rep
export NDM_SRCDIR=${NDM_ROODIR}/NDM                             # sources directory
export NDM_BUIDIR=${NDM_ROODIR}/ndm_build            # build directory

#preset_file=${NDM_SRCDIR}/cmake_files/ndm_preset_oneapi_parallel.cmake  # defines cache file with preset options
#preset_file=${NDM_SRCDIR}/cmake_files/ndm_preset_oneapi_serial.cmake  # defines cache file with preset options
preset_file=${NDM_SRCDIR}/cmake_files/ndm_lammps_oneapi.cmake  # defines cache file with preset options
#preset_file=${NDM_SRCDIR}/cmake_files/ndm_preset_gnu_parallel.cmake  # defines cache file with preset options
#preset_file=${NDM_SRCDIR}/cmake_files/ndm_preset_gnu_serial.cmake  # defines cache file with preset options
n_proc=8    # number of proc to use for compilation
























# ---------------------------- additionals options, in case cmake doesn't find packages automatically -----------------------
cmake_options=
# Lammps path and flags
#lammps_path=/home/catB/jd270899/.local/lib64/liblammps.so
#lammps_cflags=$(pkg-config liblammps --cflags)
#lammps_cflags=$(echo "$lammps_cflags" | sed 's/ //g') # remove space from the pkg-config result
#lammps_libs=$(pkg-config liblammps --libs)
#lammps_libs=$(echo "$lammps_libs" | sed 's/ //g') # remove space from the pkg-config result
#lammps_cmake_flags=" -DLAMMPS_PATH:STRING=$lammps_path" #-DLAMMPS_FLAGS:STRING=\"$lammps_cflags$lammps_libs\" "
#cmake_options+=$lammps_cmake_flags

# **************************************************************************************************************************
#   Define bash script options
# **************************************************************************************************************************

build_project=1             # build project by default
clean_flag='--clean-first'  # clean project by default
verbose_flag=''             # don't print additionnal compilation information by default
test_verbose=''             # tests build is not silent by default.
run_tests=0                 # don't run tests by default
while getopts ":h|:e|:n|:t|:v" option; do
    case $option in
        h)      # Displays this help
            echo "This script compiles and tests NDM using cmake."
            echo
            echo "Syntax : compile_ndm.sh [-h|-n|-v]"
            echo "options:"
            echo "-h    Print this help."
            echo "-n    Don't clean build directory everytime."
            echo "-t    Build and run tests."
            echo "-e    Set environment but don't run cmake and don't compile. Use with -t to only compile tests."
            echo "-v    Verbose for cmake."
            exit;;
        e)  # set environment only, use with -e option
            build_project=0 # set environment but don't run cmake
            ;;
        n)  # don't clean build directory everytime, use with -n option
            clean_flag=''     # empty to stop cleaning
            ;;
        t)  # build and run tests, use with -t option
            run_tests=1
            ;;
        v) # verbose
            verbose_flag='-v'
            ;;
    esac
done

# **************************************************************************************************************************
#   Compile NDM
# **************************************************************************************************************************
if [ $build_project -eq 1 ] # skip with -e option
then
    mkdir ${NDM_BUIDIR} 
    cd ${NDM_BUIDIR}
    rm CMakeCache.txt       #instead of --fresh (cmake < 3.24)
    cmake -C $preset_file -B ${NDM_BUIDIR} -S ${NDM_SRCDIR} 
    cmake --build ${NDM_BUIDIR} --parallel $n_proc $clean_flag $verbose_flag # build project
    clean_flag=''                   # no need to clean again for tests
    test_verbose='2>&1 > /dev/null' # make tests build silent
fi

if [ $run_tests -eq 1 ]     # use with option -t
then
    mkdir ${NDM_BUIDIR} 
    cd ${NDM_BUIDIR}
    cmake -C $preset_file -B ${NDM_BUIDIR} -S ${NDM_SRCDIR} $test_verbose
    cmake --build . --target build_tests --parallel $n_proc $clean_flag $verbose_flag              # build tests
    ctest
fi
