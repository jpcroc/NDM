#!/bin/bash -f
# **************************************************************************************************************************
#   Set environment
# **************************************************************************************************************************
# --------------------------------------- necessary environment to run or build ndm ----------------------------------------
# source ~/.bashrc
module load gcc
source /home/catB/jd270899/softwares/oneAPI/oneapi/setvars.sh --config=~/softwares/oneAPI/oneapi/config.txt
export PKG_CONFIG_PATH=${PKG_CONFIG_PATH}:/home/catB/jd270899/.local/lib64/pkgconfig
#export PKG_CONFIG_PATH=${PKG_CONFIG_PATH}:/home/catB/jd270899/.local/lib/pkgconfig

###endofsetenv - keep this line to generate setenv.sh automatically from this script
# ---------------------------------------- more environment, only relevant for build ---------------------------------------
export NDM_ROODIR=/home/catB/jd270899/Documents/NDM
#export NDM_ROODIR=/home/croc/NDM/ndm2021_cv
export NDM_SRCDIR=${NDM_ROODIR}/NDM                             # sources directory
export NDM_BUIDIR=${NDM_ROODIR}/ndm_build_ifort_para            # build directory

preset_file=${NDM_SRCDIR}/cmake_files/ndm_preset_oneapi_parallel.cmake  # defines cache file with preset options
n_proc=4    # number of proc to use for compilation

#sed '/^###endofsetenv/Q' $0 > ${NDM_BUIDIR}/setenv.sh # keeps lines above ###endofsetenv to create setenv.sh
# ---------------------------- additional options, in case cmake doesn't find packages automatically -----------------------






cmake_options=
# Lammps path and flags
#lammps_path=/home/catB/jd270899/.local/lib64/liblammps.so
#lammps_cflags=$(pkg-config liblammps --cflags)
#lammps_cflags=$(echo "$lammps_cflags" | sed 's/ //g') # remove space from the pkg-config result
#lammps_libs=$(pkg-config liblammps --libs)
#lammps_libs=$(echo "$lammps_libs" | sed 's/ //g') # remove space from the pkg-config result
#lammps_cmake_flags=" -DLAMMPS_PATH:STRING=$lammps_path" #-DLAMMPS_FLAGS:STRING=\"$lammps_cflags$lammps_libs\" "
#cmake_options+=$lammps_cmake_flags

# MKL
export MKL_ROOT=/home/catB/jd270899/softwares/oneAPI/oneapi/mkl/2023.2.0                # MKL location, use when cmake doesn't find MKL
export MKL_DIR=/home/catB/jd270899/softwares/oneAPI/oneapi/mkl/2023.2.0/lib/cmake/mkl   # MKL config location, use when cmake doesn't find MKLConfig.cmake
# cmake_options+=" -DMKL_THREADING:STRING=sequential"   # configure MKL to use relevant libraries,
# cmake_options+=" -DMKL_INTERFACE:STRING=lp64"         # find more information on available options in MKLConfig.cmake,
                                                                        # ! should be preset above in cmake cache file.
#export MPI_HOME=/home/catB/jd270899/softwares/openmpi-4.1.6/bin         # MPI location, use when cmake doesn't find MPI
# export FC=mpiifort      # set Fortran compiler, ignored if already preset in the cmake cache file !
















# **************************************************************************************************************************
#   Define bash script options
# **************************************************************************************************************************

build_project=1             # build project by default
clean_flag='--clean-first'  # clean project by default
verbose_flag=''             # don't print additionnal compilation information by default
run_tests=0                 # don't run tests by default
while getopts ":h|:e|:n|:t|:v" option; do
    case $option in
        h)      # Displays this help
            echo "This script compiles and tests NDM for gatsby using cmake."
            echo
            echo "Syntax : jinstall.para.gatsby.bash [-h|-n|-v]"
            echo "options:"
            echo "-h    Print this help."
            exit;;
        e)  # set environment only, use with -e option
            build_project=0 # set environment but don't run cmake
            ;;
        n)  # don't clean build directory everytime, use with -n option
            clean_flag=''     # empty to stop cleaning
            ;;
        t)  # build and run tests only, use with -t option
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
    cd ${NDM_BUIDIR}   # important
    cmake $cmake_options -C $preset_file -B ${NDM_BUIDIR} -S ${NDM_SRCDIR} --fresh                 # generate project
    cmake --build ${NDM_BUIDIR} --parallel $n_proc $clean_flag $verbose_flag # build project
fi

if [ $run_tests -eq 1 ]     # use with option -t
then
    cmake --build . --target test               # run tests
    # # # make an integration example...
    rm -rf  ./tmp
    cp -rf ${NDM_SRCDIR}/examples/ndm_mpi ./tmp
    cd ./tmp
#    mpirun -n $n_proc ${NDM_BUIDIR}/bin/ndm_main.exe
fi
