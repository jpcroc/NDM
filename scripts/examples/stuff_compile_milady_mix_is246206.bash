#!/usr/bin/env bash

# on is246206 compile milady and run ctest and run python tests
# rename and modification and use at your own risk

# bash
# export PATH=/export/home/catA/wambeke/MLD/MILADY/scripts:${PATH}
# source /export/home/catA/wambeke/MLD/MILADY/scripts/examples/stuff_compile_milady_mix_is246206.bash
# source /export/home/catA/wambeke/MLD/MILADY/scripts/examples/stuff_compile_milady_mix_is246206.bash | grep -e '-- ..' | grep MLD_MKL_LIB
# f_stuff_mix_is246206

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

function f_stuff_mix_is246206 {

  # to stop stuff sometimes
  envs --create_setenv="${HOME}/mld_setenv_mix.tmp" # get current bash environ status
  export mld_do_python_tests=ON   # to do python integration tests
  export mld_do_make=ON           # to compile milady

  export MLD_OPT_TRAC=OFF          # to log cmake all variables contents

  # check bash environment is clean and previously coherent
  [ -z ${MLD_TYPE} ] && export MLD_TYPE=MIX  # mandatory else risky
  [ ! "${MLD_TYPE}" == "MIX" ] && in_red_pause "ERROR: bash environment is not clean, continue at your risks"

  # memorize origin clean env vars
  [ -z ${LD_LIBRARY_PATH_ORIG} ] && export LD_LIBRARY_PATH_ORIG=${LD_LIBRARY_PATH}
  [ -z ${PATH_ORIG} ] && export PATH_ORIG=${PATH}

  # fixed may be as 'alias cmake=cmake3'
  tmp=$(cmake --version | head -1)  # have to be cmake3 version > 3.x.y (3.17.5)
  [[ ${tmp} = *" 3."* ]] || in_red_pause 'ERROR: '${tmp}' , 3.x.y is mandatory'
  in_green "using cmake3"

  # for example. out-of-source-directory build and install
  # policy should be: stay on these useful environment names !(keep names, change values)
  in_green "f_stuff set MLD_ etc."
  export LOC_ROODIR="/export/home/catA/wambeke"
  export MLD_OAP_ROODIR="/export/home/catA/intel"
  export MLD_OAP_INSDIR="${MLD_OAP_ROODIR}/oneapi"
  export MLD_ROODIR="${LOC_ROODIR}/MLD"
  export MLD_SRCDIR="${MLD_ROODIR}/MILADY"
  export MLD_BUIDIR="${MLD_ROODIR}/mld_build"     # is user choice
  export MLD_INSDIR="${MLD_ROODIR}/mld_install"   # is user choice
  export MLD_TESDIR="${MLD_ROODIR}/mld_tests"     # is user choice
  export MLD_MPI_ROODIR="${MLD_ROODIR}/mpi"  # /usr/local/iopenmpi.4.1.2/lib64/
  export MLD_MPI_INSDIR="${MLD_MPI_ROODIR}/install_412_intel"  # /usr/local/iopenmpi.4.1.2/lib64/
  export MLD_CMAKE_BUILD_TYPE=${MLD_CMAKE_BUILD_TYPE:-Release}  # or "Debug" #
  export MLD_SETENV_DONE=ON

  which envs > /dev/null || export PATH=${MLD_SRCDIR}/scripts:${PATH}
  # in_green 'which envs as '$(which envs)
  envs -p MLD_
  envs -p mld_

  [ ! -d ${MLD_OAP_ROODIR} ] && in_red_pause 'problem: install oneAPI as ${MLD_OAP_ROODIR}' ${MLD_OAP_ROODIR}
  [ ! -d ${MLD_OAP_INSDIR} ] && in_red_pause 'problem: install oneAPI as ${MLD_OAP_ROODIR}/oneapi' ${MLD_OAP_INSDIR}
  [ ! -d ${MLD_ROODIR} ] && in_red_pause 'problem: inexisting ${MLD_ROODIR}='${MLD_ROODIR}'.\ntry: mkdir ${MLD_ROODIR}'
  [ ! -d ${MLD_SRCDIR} ] && in_red_pause 'problem: inexisting ${MLD_SRCDIR}='${MLD_SRCDIR}'.\ntry: cd ${MLD_ROODIR}; git clone https://codev-tuleap.intra.cea.fr/plugins/git/milady/MiLaDy.git MILADY'
  [ ! -d ${MLD_TESDIR} ] && in_red_pause 'problem: inexisting ${MLD_TESDIR}='${MLD_TESDIR}'\nfix it.'

  alias cd_loc="cd ${LOC_ROODIR} ; pwd"
  alias cd_milady="cd ${MLD_SRCDIR} ; pwd"
  alias cd_mld="cd ${MLD_ROODIR} ; pwd"
  alias cd_tests="cd ${MLD_TESDIR} ; pwd"

  # in_green "*** set oneapi as default choice ***"
  # source ${MLD_OAP_ROODIR}/setvars.sh
  # export MLD_OAP_ROODIR=.../intel
  # export MLD_OAP_INSDIR=${MLD_OAP_ROODIR}/oneapi        # oneAPI install
  # export MLD_ROODIR=$(pwd)  # for example...

  # case Cosmin linux suze have not module command
  # source ${MLD_OAP_INSDIR}/setvars.sh  # but on config file

  if [ -z ${SETVARS_COMPLETED} ]; then
    # Intel Trace Analyzer and Collector (ITAC) is an MPI profiling and tracing tool
    in_red 'need launched setvars.sh oneAPI, set it for MIX. compiler and mkl and itac only'
    export CONFIG_SETVARS=${MLD_ROODIR}/config_setvars_oneapi.tmp
    cat <<EOT > ${CONFIG_SETVARS}
default=exclude
compiler=latest
mkl=latest
itac=latest
EOT
    in_green 'config oneAPI setvars on ${MLD_ROODIR}/config_setvars_oneapi.tmp'
    cat ${CONFIG_SETVARS} # to verify
    # if may be done yet:
    # WARNING: setvars.sh has already been run. Skipping re-execution.
    # To force a re-execution of setvars.sh, use the '--force' option.
    # Using '--force' can result in excessive use of your environment variables.
    source ${MLD_OAP_INSDIR}/setvars.sh --config=${CONFIG_SETVARS}

    # suze have not module: this is obsolete
    # module use ${MLD_OAP_INSDIR}/modulefiles
    # module avail
    # module load compiler/latest
    # module load mkl/latest
  else
    in_red 'setvars.sh and ${SETVARS_COMPLETED}' ${SETVARS_COMPLETED}
    in_red 'launched setvars.sh oneAPI done yet, hope it is OK (without intel mpi)'
    # hope there is NOT intel mpi
    envs -g mpi
  fi

  # some verifications
  [ -z ${MKLROOT} ] && in_red_pause 'problem: inexisting ${MKLROOT} from oneAPI\nfix it.'

  # set by oneAPI setvars.sh
  envs -e MKLROOT
  envs -g oneapi/itac
  envs -g oneapi/compiler

  tmp=$(which ifort 2> /dev/null) && in_green 'ifort at '${tmp} || in_red_pause 'problem no reached ifort compiler'


  cd ${MLD_SRCDIR}/scripts
  source ./compile_milady.bash

  # envs -e PATH,LD_LIBRARY_PATH

  unset MLD_NPROC # let it automatic
  # export MLD_NPROC=10
  # export MLD_OPT_TRAC=ON
  f_compile_milady_mix # do cmake

  if [ "${mld_do_make}" == "OFF" ]; then  # stop here ?
    cd ${MLD_ROODIR}
    return
  fi
  f_make # do make
  in_green "make install"
  make install

  in_green 'set ${MILADYPY_ROOT_DIR} on miladypy'
  export MILADYPY_ROOT_DIR=${MLD_SRCDIR}
  in_green 'set PATH & LD_LIBRARY_PATH on openmpi'
  export PATH=${MLD_MPI_INSDIR}/bin:${PATH}
  export LD_LIBRARY_PATH=${MLD_MPI_INSDIR}/lib:${LD_LIBRARY_PATH}

  #if [ ! -z ${mld_do_obsolete} ]; then
    # created a file mld_setenv.bash from cmake to resume/fix next
    #in_green 'set PATH & LD_LIBRARY_PATH on mkl itac compiler'
    #export LD_LIBRARY_PATH=${MKL_ROOT}/lib/intel64:${LD_LIBRARY_PATH}
    # export LD_LIBRARY_PATH=${MLD_OAP_INSDIR}/itac/2021.5.0/bin/rtlib:${LD_LIBRARY_PATH}
    #export LD_LIBRARY_PATH=${VT_ROOT}/bin/rtlib:${LD_LIBRARY_PATH}
    # export LD_LIBRARY_PATH=${MLD_OAP_INSDIR}/compiler/2022.0.1/linux/compiler/lib/intel64_lin:${LD_LIBRARY_PATH}
    #export LD_LIBRARY_PATH=${CMPLR_ROOT}/linux/compiler/lib/intel64_lin:${LD_LIBRARY_PATH}
  #fi

  in_green "create_diff_setenv from "${HOME}/mld_setenv.tmp
  envs --create_diff_setenv="${HOME}/mld_setenv_mix.tmp" # get current bash environ status

  envs -e LD_LIBRARY_PATH,PATH

  # example find libmpi_mpifh.so
  # in_green "find libmpi_mpifh.so"
  # locate libmpi_mpifh.so
  # ldconfig -p | grep libmpi_mpifh.so

  # f_explore_milady   # long listing only for debug
  f_ctest_milady

  if [ "${mld_do_python_tests}" == "OFF" ]; then  # stop here ?
    cd ${MLD_ROODIR}
    return
  fi

  in_green "set python3 for python small_tests"
  f_conda_init  # function to set environment conda/python at is246206
  conda env list
  conda activate py3qt5

  # options for mld unittest python
  mld_plotFailed=ON
  mld_verboseTest=ON
  f_python_tests_milady
}


in_green 'to compile milady and execute tests, type:\n  f_stuff_mix_is246206'
[ "${USER}" == "wambeke" ] && f_stuff_mix_is246206   # do it for me
