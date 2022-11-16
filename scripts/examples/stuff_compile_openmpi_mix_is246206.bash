#!/usr/bin/env bash

# on is246206 compile openmpi milady
# rename and modification and use at your own risk

# bash
# export PATH=/export/home/catA/wambeke/MLD/MILADY/scripts:${PATH}
# source /export/home/catA/wambeke/MLD/MILADY/scripts/examples/stuff_compile_openmpi_mix_is246206.bash
# f_stuff_is246206_mix

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

function f_stuff_openmpi_is246206_mix {

  # to get a setenv diff file
  envs --create_setenv="${HOME}/mld_setenv_openmpi.tmp" # get current bash environ status

  # check bash environment is clean and previously coherent
  [ -z ${MLD_TYPE} ] && export MLD_TYPE=MIX  # mandatory else risky
  [ ! "${MLD_TYPE}" == "MIX" ] && in_red_pause "ERROR: bash environment is not clean, continue at your risks"

  # memorize origin clean env vars
  [ -z ${LD_LIBRARY_PATH_ORIG} ] && export LD_LIBRARY_PATH_ORIG=${LD_LIBRARY_PATH}
  [ -z ${PATH_ORIG} ] && export PATH_ORIG=${PATH}

  # for example. out-of-source-directory build and install
  # policy should be: stay on these useful environment names !(keep names, change values)
  in_green "f_stuff set MLD_ etc."
  export LOC_ROODIR="/export/home/catA/wambeke"
  export MLD_OAP_ROODIR="/export/home/catA/intel"
  export MLD_OAP_INSDIR="${MLD_OAP_ROODIR}/oneapi"
  export MLD_ROODIR="${LOC_ROODIR}/MLD"
  export MLD_SRCDIR="${MLD_ROODIR}/MILADY"
  #export MLD_BUIDIR="${MLD_ROODIR}/mld_build"     # is user choice
  #export MLD_INSDIR="${MLD_ROODIR}/mld_install"   # is user choice
  #export MLD_TESDIR="${MLD_ROODIR}/mld_tests"     # is user choice
  export MLD_MPI_INSDIR="${MLD_ROODIR}/mpi/install_412_intel"  # /usr/local/iopenmpi.4.1.2/lib64/
  #export MLD_CMAKE_BUILD_TYPE=${MLD_CMAKE_BUILD_TYPE:-Release}  # or "Debug" #
  export MLD_SETENV_DONE=ON
  export MLD_NPROC=${MLD_NPROC:-8}

  which envs > /dev/null || export PATH=${MLD_SRCDIR}/scripts:${PATH}
  # in_green 'which envs as '$(which envs)
  envs -p MLD_
  # envs -p mld_

  [ ! -d ${MLD_OAP_ROODIR} ] && in_red_pause 'problem: install oneAPI as ${MLD_OAP_ROODIR}' ${MLD_OAP_ROODIR}
  [ ! -d ${MLD_OAP_INSDIR} ] && in_red_pause 'problem: install oneAPI as ${MLD_OAP_ROODIR}/oneapi' ${MLD_OAP_INSDIR}
  [ ! -d ${MLD_ROODIR} ] && in_red_pause 'problem: inexisting ${MLD_ROODIR}='${MLD_ROODIR}'.\ntry: mkdir ${MLD_ROODIR}'
  #[ ! -d ${MLD_SRCDIR} ] && in_red_pause 'problem: inexisting ${MLD_SRCDIR}='${MLD_SRCDIR}'.\ntry: cd ${MLD_ROODIR}; git clone https://codev-tuleap.intra.cea.fr/plugins/git/milady/MiLaDy.git MILADY'
  #[ ! -d ${MLD_TESDIR} ] && in_red_pause 'problem: inexisting ${MLD_TESDIR}='${MLD_TESDIR}'\nfix it.'

  alias cd_loc="cd ${LOC_ROODIR} ; pwd"
  #alias cd_milady="cd ${MLD_SRCDIR} ; pwd"
  alias cd_mld="cd ${MLD_ROODIR} ; pwd"
  #alias cd_tests="cd ${MLD_TESDIR} ; pwd"

  # case Cosmin linux suze have not module command
  # source ${MLD_OAP_INSDIR}/setvars.sh  # but on config file

  if [ -z ${SETVARS_COMPLETED} ]; then
    # Intel Trace Analyzer and Collector (ITAC) is an MPI profiling and tracing tool
    in_red 'need launched setvars.sh oneAPI, set it for MIX. compiler only (and no mkl itac)'
    export CONFIG_SETVARS=${MLD_ROODIR}/config_setvars_oneapi.tmp
    cat <<EOT > ${CONFIG_SETVARS}
default=exclude
compiler=latest
EOT
    in_green 'config oneAPI setvars on ${MLD_ROODIR}/config_setvars_oneapi.tmp'
    cat ${CONFIG_SETVARS} # to verify
    # if may be done yet:
    # WARNING: setvars.sh has already been run. Skipping re-execution.
    # To force a re-execution of setvars.sh, use the '--force' option.
    # Using '--force' can result in excessive use of your environment variables.
    source ${MLD_OAP_INSDIR}/setvars.sh --config=${CONFIG_SETVARS}
  else
    in_red 'setvars.sh and ${SETVARS_COMPLETED}' ${SETVARS_COMPLETED}
    in_red 'launched setvars.sh oneAPI done yet, hope it is OK (without intel mpi)'
    # hope there is NOT intel mpi
    envs -g '/mpi' || n_red_pause 'problem: needed there is NOT intel mpi'
  fi

  # some verifications
  #hope there is NOT intel mpi
  envs -g '/mpi' || n_red_pause 'problem: needed there is NOT intel mpi'
  # [ -z ${MKLROOT} ] && in_red_pause 'problem: inexisting ${MKLROOT} from oneAPI\nfix it.'

  # set by oneAPI setvars.sh
  #envs -e MKLROOT
  #envs -g oneapi/itac
  envs -g oneapi/compiler

  tmp=$(which icc 2> /dev/null) && in_green 'icc at '${tmp} || in_red_pause 'problem no reached icc compiler'
  tmp=$(which icpc 2> /dev/null) && in_green 'icpc at '${tmp} || in_red_pause 'problem no reached icpc compiler'
  tmp=$(which ifort 2> /dev/null) && in_green 'ifort at '${tmp} || in_red_pause 'problem no reached ifort compiler'

  cd ${MLD_SRCDIR}/scripts
  source ./compile_openmpi.bash

  envs -e PATH,LD_LIBRARY_PATH

  unset MLD_NPROC # let it automatic
  # export MLD_NPROC=10
  f_setenv_openmpi_412_wambeke
  f_compile_openmpi_intel  # compile openmpi with intel compilers ... for milady mix

  export PATH=${MLD_MPI_INSDIR}/bin:${PATH}
  export LD_LIBRARY_PATH=${MLD_MPI_INSDIR}/lib:${LD_LIBRARY_PATH}
  envs --create_diff_setenv="${HOME}/mld_setenv_openmpi.tmp" # get current bash environ status

  tree -d ${MLD_MPI_INSDIR} # check
  [ -f "${MLD_MPI_INSDIR}/bin/mpirun" ] \
    && in_green 'OK f_stuff_openmpi_is246206_mix compile openmpi with intel compilers' \
    || in_red_pause 'KO f_stuff_openmpi_is246206_mix compile openmpi with intel compilers'
}


in_green 'to compile openmpi, type:\n  f_stuff_openmpi_is246206_mix'
[ "${USER}" == "wambeke" ] && f_stuff_openmpi_is246206_mix   # do it for me
