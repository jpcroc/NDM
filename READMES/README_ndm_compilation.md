
## NDM GNU serial or parallel compilation on local user host linux


### example set opempi

Do it Before

```
module avail
  ---------- /usr/share/Modules/modulefiles
  dot  module-git  module-info  modules  null  use.own  
  ---------- /usr/share/modulefiles
  mp-x86_64  mpi/mpich-x86_64  mpi/openmpi-x86_64

bash
module load mpi/openmpi-x86_64
```
.

### example compile GNU serial or parallel

```
hostname  # is246206.intra.cea.fr
whoami    # wambeke
export NDM_ROODIR=/export/home/catA/wambeke/NDM_ROODIR
cd ${NDM_ROODIR}
# git clone --branch ndm2021_cv ssh://gitolite@ssh-codev-tuleap.intra.cea.fr:2044/ndm/NDM.git NDM

export NDM_SRCDIR=${NDM_ROODIR}/NDM
export NDM_BUIDIR=${NDM_ROODIR}/ndm_build_gnu
export NDM_INSDIR=${NDM_ROODIR}/ndm_install_gnu

bash   # useful to rewind clean environment with exit

export PATH=${NDM_SRCDIR}/scripts:${PATH}  # to get NMD/scripts
envs -p NDM_

# preset.cmakes files makes this important fortran compiler choice useless
unset FC

cd ${NDM_ROODIR}
rm -rf ${NDM_BUIDIR}
rm -rf ${NDM_INSDIR}
mkdir ${NDM_BUIDIR}

cd ${NDM_BUIDIR}   # important

###### GNU serial
cmake -C ${NDM_SRCDIR}/cmake_files/ndm_preset_gnu_serial.cmake -S ${NDM_SRCDIR}

###### GNU parallel need openmpi, do that one time only
if [ -z "${MPI_BIN}" ] ; then
  module avail
  module load mpi/openmpi-x86_64
fi

cmake -C ${NDM_SRCDIR}/cmake_files/ndm_preset_gnu_parallel.cmake -S ${NDM_SRCDIR}

ccmake .  # use it only for display

# nproc     # --> 20 is246206
# mate-system-monitor &

cmake --build . --target all --parallel 15
cmake --build . --target install
cmake --build . --target test

# checks
tree -d ${NDM_BUIDIR} # ../ndm_build_gnu/
tree -d ${NDM_INSDIR} # ../ndm_install_gnu/
# tree -d ../ndm_*gnu
ldd ${NDM_BUIDIR}/bin/*

# make an integration example...
rm -rf  ./tmp
cp -rf ${NDM_SRCDIR}/examples/ndm_lammps_serial ./tmp
cd ./tmp

# ndm_preset_gnu_serial
${NDM_INSDIR}/bin/ndm_main.exe
---> ok

# parallel
mpirun -n 4 ${NDM_INSDIR}/bin/ndm_main.exe
---> ok

```
.

## NDM mpi intel oneApi serial or parallel compilation on local user host linux

.

### example set oneApi serial or parallel

Do it Before

- useful to launch setvars.sh oneAPI.
- set it for INTEL compilers and mkl and itac only.
- Intel Trace Analyzer and Collector (ITAC) is an MPI profiling and tracing tool.


```
# oneapi installation dir for is246206
export NDM_OAP_ROODIR=/export/home/catA/intel
export NDM_OAP_INSDIR=${NDM_OAP_ROODIR}/oneapi

export NDM_ROODIR=/export/home/catA/wambeke/NDM_ROODIR
export NDM_SRCDIR=${NDM_ROODIR}/NDM
export NDM_BUIDIR=${NDM_ROODIR}/ndm_build_intel
export NDM_INSDIR=${NDM_ROODIR}/ndm_install_intel

cd ${NDM_ROODIR}

bash   # useful to rewind clean environment with exit

export PATH=${NDM_SRCDIR}/scripts:${PATH}  # to get NMD/scripts
envs -p NDM_

if [ -z ${SETVARS_COMPLETED} ]; then    # do it only one time

  export CONFIG_SETVARS=${NDM_ROODIR}/config_setvars_oneapi.tmp

  ###### INTEL parallel need openmpi, do that one time only
  cat <<EOT > ${CONFIG_SETVARS}
default=exclude
compiler=latest
mkl=latest
mpi=latest
itac=latest
EOT

  ###### INTEL serial do not need openmpi, do that one time only
  cat <<EOT > ${CONFIG_SETVARS}
default=exclude
compiler=latest
mkl=latest
itac=latest
EOT

  echo 'config oneAPI setvars on ${NDM_ROODIR}/config_setvars_oneapi.tmp'
  cat ${CONFIG_SETVARS} # to verify
  source ${NDM_OAP_INSDIR}/setvars.sh --config=${CONFIG_SETVARS}
  echo "SETVARS_COMPLETED==${SETVARS_COMPLETED}"
  echo "I_MPI_ROOT=${I_MPI_ROOT}"  # set
  echo "MPI_HOME=${MPI_HOME}"  # not set

else

  echo 'setvars.sh and ${SETVARS_COMPLETED}' ${SETVARS_COMPLETED}
  echo 'launched setvars.sh oneAPI done yet, hope it is OK (with or without intel mpi, mkl...)'
  # hope there is intel mpi
  envs -g mpi

fi
```
.

OS suze, for example have not `module` command, this intel documented way is obsolete:

```
module use ${NDM_OAP_INSDIR}/modulefiles
module avail
module load compiler/latest
module load mkl/latest
```
.

If done yet you get:

```
WARNING: setvars.sh has already been run. Skipping re-execution.
To force a re-execution of setvars.sh, use the '--force' option.
Using '--force' can result in excessive use of your environment variables.
```
.

### example compile INTEL


```
# preset.cmakes files makes this important fortran compiler choice useless
unset FC

cd ${NDM_ROODIR}
rm -rf ${NDM_BUIDIR}
rm -rf ${NDM_INSDIR}
mkdir ${NDM_BUIDIR}

cd ${NDM_BUIDIR}   # important

###### intel serial
cmake -C ${NDM_SRCDIR}/cmake_files/ndm_preset_oneapi_serial.cmake -S ${NDM_SRCDIR}

###### intel parallel
cmake -C ${NDM_SRCDIR}/cmake_files/ndm_preset_oneapi_parallel.cmake -S ${NDM_SRCDIR}

ccmake .  # use it only for display

# nproc     # --> 20 is246206
# mate-system-monitor &

cmake --build . --target all --parallel 15
cmake --build . --target install
cmake --build . --target test


# checks
tree -d ${NDM_BUIDIR} # ../ndm_build_gnu/
tree -d ${NDM_INSDIR} # ../ndm_install_gnu/
# tree -d ../ndm_*gnu
ldd ${NDM_BUIDIR}/bin/*

# make an integration example...
rm -rf  ./tmp
cp -rf ${NDM_SRCDIR}/examples/ndm_lammps_serial ./tmp
cd ./tmp

# serial
${NDM_INSDIR}/bin/ndm_main.exe
---> ok

# parallel
mpirun -n 4 ${NDM_INSDIR}/bin/ndm_main.exe
---> ok

```
.
