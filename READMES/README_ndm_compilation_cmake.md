
## NDM compilation serial (no mpi) on on local user host linux

```
hostname  # is246206.intra.cea.fr
whoami    # wambeke
export NDM_ROODIR=/export/home/catA/wambeke/NDM_ROODIR
cd ${NDM_ROODIR}
git clone --branch ndm2021_cv ssh://gitolite@ssh-codev-tuleap.intra.cea.fr:2044/ndm/NDM.git NDM

export NDM_SRCDIR=${NDM_ROODIR}/NDM
export NDM_BUIDIR=${NDM_ROODIR}/ndm_build_gnu
export NDM_INSDIR=${NDM_ROODIR}/ndm_install_gnu


bash
export PATH=${NDM_SRCDIR}/scripts:${PATH}  # to get NMD/scripts
envs -p NDM_

cd ${NDM_ROODIR} ; rm -rf ${NDM_BUIDIR} ; rm -rf ${NDM_INSDIR} ; mkdir ${NDM_BUIDIR} ; cd ${NDM_BUIDIR}

cmake -C ${NDM_SRCDIR}/cmake_files/ndm_preset_gnu_serial.cmake ${NDM_SRCDIR}

cmake-gui ../cmake

# --parallel 4
# nproc --> 20
# mate-system-monitor &

cmake --build . --target all --parallel 15
cmake --build . --target install
cmake --build . --target test

tree -d ${NDM_BUIDIR} # ../ndm_build_gnu/
tree -d ${NDM_INSDIR} # ../ndm_install_gnu/
tree -d ../ndm_*gnu
ldd ${NDM_INSDIR}/bin/*

# make an integration example...
rm -rf  ./tmp
cp -rf ${NDM_SRCDIR}/examples/ndm_lammps_serial ./tmp
cd ./tmp
${NDM_INSDIR}/bin/ndm_main.exe
---> ok

```


### for future with opempi

TODO

```
module avail
  ---------- /usr/share/Modules/modulefiles
  dot  module-git  module-info  modules  null  use.own  
  ---------- /usr/share/modulefiles
  mp-x86_64  mpi/mpich-x86_64  mpi/openmpi-x86_64

# for gnu is useful to get mpi
bash
module load mpi/openmpi-x86_64

cmake -C ${NDM_SRCDIR}/cmake_files/ndm_preset_gnu_serial.cmake \
  -D CMAKE_VERBOSE_MAKEFILE=off \
  -D CMAKE_INSTALL_PREFIX=${NDM_INSDIR} \
  -D ENABLE_TESTING=on \
  ${NDM_SRCDIR}

```
