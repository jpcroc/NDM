#!/usr/bin/env bash

# on is246206 compile openmpi milady
# rename and modification and use at your own risk

# bash
# export PATH=/export/home/catA/wambeke/MLD/MILADY/scripts:${PATH}
# source /export/home/catA/wambeke/MLD/MILADY/scripts/examples/stuff_compile_openmpi_gnu_is246206.bash
# f_stuff_gnu_is246206

function in_red () {
  # Red bold
  echo -e "\e[31m\e[1m"$@"\e[0m"
}

function in_red_pause () {
  # Red bold
  echo -e "\e[31m\e[1m"$@"\e[0m"
  f_pause
}

function f_pause {
  read -p 'Continue? (ctrl-c to stop)' rep
}

function in_green () {
  # green bold
  echo -e "\e[32m\e[1m"$@"\e[0m"
}

function f_stuff_openmpi_gnu_is246206 {

  # to get a setenv diff file
  envs --create_setenv="${HOME}/mld_setenv_openmpi.tmp" # get current bash environ status

  # check bash environment is clean and previously coherent
  [ -z ${MLD_TYPE} ] && export MLD_TYPE=GNU  # mandatory else risky
  [ ! "${MLD_TYPE}" == "GNU" ] && in_red_pause "ERROR: bash environment is not clean, continue at your risks"

  # memorize origin clean env vars
  [ -z ${LD_LIBRARY_PATH_ORIG} ] && export LD_LIBRARY_PATH_ORIG=${LD_LIBRARY_PATH}
  [ -z ${PATH_ORIG} ] && export PATH_ORIG=${PATH}

  # for example. out-of-source-directory build and install
  # policy should be: stay on these useful environment names !(keep names, change values)
  in_green "f_stuff set MLD_ etc."
  export LOC_ROODIR="/export/home/catA/wambeke"
  export MLD_ROODIR="${LOC_ROODIR}/MLD"
  export MLD_SRCDIR="${MLD_ROODIR}/MILADY"
  export MLD_MPI_ROODIR="${MLD_ROODIR}/mpi"
  export MLD_MPI_INSDIR="${MLD_MPI_ROODIR}/install_412_intel"  # /usr/local/iopenmpi.4.1.2/lib64/
  export MLD_SETENV_DONE=ON
  export MLD_NPROC=${MLD_NPROC:-8}

  which envs > /dev/null || export PATH=${MLD_SRCDIR}/scripts:${PATH}
  # in_green 'which envs as '$(which envs)
  envs -p MLD_
  # envs -p mld_

  [ ! -d ${MLD_ROODIR} ] && in_red_pause 'problem: inexisting ${MLD_ROODIR}='${MLD_ROODIR}'.\ntry: mkdir ${MLD_ROODIR}'

  alias cd_loc="cd ${LOC_ROODIR} ; pwd"
  alias cd_mld="cd ${MLD_ROODIR} ; pwd"

  if [ ! -z ${SETVARS_COMPLETED} ]; then
    in_red 'NO need launched setvars.sh oneAPI, start from clean bash as GNU compiler only (and NO oneAPI intel)'
    return
  else
    in_green 'GNU compiler use, hope it is OK (without intel mpi)'
  fi

  # some verifications
  # hope there is NOT intel mpi
  envs -g '/mpi' || n_red_pause 'problem: needed there is NOT intel mpi'

  tmp=$(which cc 2> /dev/null) && in_green 'cc at '${tmp} || in_red_pause 'problem no reached cc compiler'
  tmp=$(which gcc 2> /dev/null) && in_green 'gcc at '${tmp} || in_red_pause 'problem no reached gcc compiler'
  tmp=$(which gfortran 2> /dev/null) && in_green 'gfortran at '${tmp} || in_red_pause 'problem no reached gfortran compiler'


  cd ${MLD_SRCDIR}/scripts
  source ./compile_openmpi.bash

  envs -e PATH,LD_LIBRARY_PATH

  # unset MLD_NPROC # let it automatic
  # export MLD_NPROC=10
  f_setenv_openmpi_412_wambeke
  f_compile_openmpi_gnu  # compile openmpi with intel compilers ... for milady gnu

  export PATH=${MLD_MPI_INSDIR}/bin:${PATH}
  export LD_LIBRARY_PATH=${MLD_MPI_INSDIR}/lib:${LD_LIBRARY_PATH}
  envs --create_diff_setenv="${HOME}/mld_setenv_openmpi.tmp" # get current bash environ status

  tree -d ${MLD_MPI_INSDIR} # check
  [ -f "${MLD_MPI_INSDIR}/bin/mpirun" ] \
    && in_green 'OK f_stuff_openmpi_gnu_is246206 compile openmpi with intel compilers' \
    || in_red_pause 'KO f_stuff_openmpi_gnu_is246206 compile openmpi with intel compilers'
}


in_green 'to compile openmpi, type:\n  f_stuff_openmpi_gnu_is246206'
[ "${USER}" == "wambeke" ] && f_stuff_openmpi_gnu_is246206   # do it for me
