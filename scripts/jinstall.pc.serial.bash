#!/bin/bash -f
export NDM_ROODIR=/mnt/c/Users/jc148490/lin/DM/codesndm/CVW3
cd ${NDM_ROODIR}

export NDM_SRCDIR=${NDM_ROODIR}/NDM
export NDM_BUIDIR=${NDM_ROODIR}/ndm_build_gnu
export NDM_INSDIR=${NDM_ROODIR}/ndm_install_gnu

export PATH=${NDM_SRCDIR}/scripts:${PATH}  # to get NMD/scripts
envs -p NDM_

function f_cmake_clean_doc {
  rm -rf ${NDM_SRCDIR}/docs/doc_2023/_build
  rm -rf ${NDM_SRCDIR}/docs/doc_2023/_tmp
  rm -rf ${NDM_SRCDIR}/docs/doc_2023/__pycache__
  tree -d ${NDM_SRCDIR}/docs
}

function f_cmake_clean {
  f_cmake_clean_doc
  cd ${NDM_ROODIR}
  rm -rf ${NDM_BUIDIR}
  rm -rf ${NDM_INSDIR}
  mkdir ${NDM_BUIDIR}
  cd ${NDM_BUIDIR}   # important
}

# preset.cmakes files makes this important fortran compiler choice useless
unset FC

f_cmake_clean   # important

# you could append 'cmake -DCMAKE_VERBOSE_MAKEFILE=on ...''

# example compile with trace, only 'jupy'
# cmake -D NDM_OPT_TRACE=ON -C ${NDM_SRCDIR}/cmake_files/ndm_preset_gnu_serial.cmake -S ${NDM_SRCDIR} 2>&1 | grep jupy

# example compile doc_2023
# cmake -D NDM_OPT_COMPILE_DOC=ON -C ${NDM_SRCDIR}/cmake_files/ndm_preset_gnu_serial.cmake -S ${NDM_SRCDIR} ${NDM_SRCDIR}

# compile standart (without doc and trace)

###### GNU serial
cmake -C ${NDM_SRCDIR}/cmake_files/ndm_preset_gnu_serial.cmake -S ${NDM_SRCDIR}

###### GNU parallel need openmpi, do that one time only
#PROBLEME


  # C++ = mpic++
  # CC  = mpicc
  # F77 = mpif77  -cpp # -warn all -check all
  # F90 = mpif90  -cpp # -warn all -check all


  # MPI_LIB  = -lmpi
# MKL_LIB  = -L/mnt/c/Users/jc148490/lin/INFO/lapack/lapack-3.10.0 -llapack -lrefblas

# if [ -z "${MPI_BIN}" ] ; then   C   MPI_BIN ???
  # module avail
  # module load mpi/openmpi-x86_64
# fi
 #cmake -C ${NDM_SRCDIR}/cmake_files/ndm_preset_gnu_parallel.cmake -S ${NDM_SRCDIR}

# ccmake .  # use it only for display

# nproc     # --> 20 is246206
# mate-system-monitor &

cmake --build . --target all --parallel 4
cmake --build . --target install
cmake --build . --target test

# checks
tree -d ${NDM_BUIDIR} # ../ndm_build_gnu/
ldd ${NDM_BUIDIR}/bin/*

# # # make an integration example...
 rm -rf  ./tmp
 cp -rf ${NDM_SRCDIR}/examples/ndm_lammps_serial ./tmp
 cd ./tmp

${NDM_INSDIR}/bin/ndm_main.exe
# ---> ok

