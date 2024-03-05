#!/bin/bash -f
export NDM_ROODIR=/ccc/cont002/home/den/croc/CVW3B/CVW3
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

###### GNU serial
cmake -C ${NDM_SRCDIR}/cmake_files/ndm_preset_gnu_serial.cmake -S ${NDM_SRCDIR}


cmake --build . --target all --parallel 4
#cmake --build . --target install
cmake --build . --target test

# checks
tree -d ${NDM_BUIDIR} # ../ndm_build_gnu/
#ldd ${NDM_BUIDIR}/bin/*

# # # make an integration example...
 rm -rf  ./tmp
 cp -rf ${NDM_SRCDIR}/examples/ndm_lammps_serial ./tmp
 cd ./tmp

#ndm_preset_gnu_serial #????#
 ${NDM_BUIDIR}/bin/ndm_main.exe
# ---> ok



