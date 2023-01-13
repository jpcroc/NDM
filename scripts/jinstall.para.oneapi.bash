#!/bin/bash -f
#source ~/.bashrc
module purge
module load inteloneapi/21.4.0
module load mpi/intelmpi/21.4.0

export NDM_ROODIR=/ccc/cont002/home/den/croc/CVW3
cd ${NDM_ROODIR}

export NDM_SRCDIR=${NDM_ROODIR}/NDM
export NDM_BUIDIR=${NDM_ROODIR}/ndm_build_ifort_para
export NDM_INSDIR=${NDM_ROODIR}/ndm_install_ifort_para

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

# which mpif90

export MPI_HOME=/ccc/products/intelmpi-21.4.0/system/default/21.4.0

f_cmake_clean   # important


cmake -C ${NDM_SRCDIR}/cmake_files/ndm_preset_oneapi_parallel.cmake -S ${NDM_SRCDIR}

cmake --build . --target all --parallel 4
#cmake --build . --target install
cmake --build . --target test


# checks
tree -d ${NDM_BUIDIR} 
#tree -d ${NDM_INSDIR} 
# tree -d ../ndm_*gnu
#ldd ${NDM_BUIDIR}/bin/*


# # # make an integration example...
 rm -rf  ./tmp
 cp -rf ${NDM_SRCDIR}/examples/ndm_mpi ./tmp
 cd ./tmp


#parallel
 mpirun -n 4 ${NDM_BUIDIR}/bin/ndm_main.exe
# ---> ok

