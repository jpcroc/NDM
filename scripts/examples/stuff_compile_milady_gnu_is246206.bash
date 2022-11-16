#!/usr/bin/env bash

# on is246206 compile milady and run ctest and run python tests
# rename and modification and use at your own risk

# bash
# export PATH=/export/home/catA/wambeke/MLD/MILADY/scripts:${PATH}
# source /export/home/catA/wambeke/MLD/MILADY/scripts/examples/stuff_compile_milady_gnu_is246206.bash
# source /export/home/catA/wambeke/MLD/MILADY/scripts/examples/stuff_compile_milady_gnu_is246206.bash | grep -e '-- ..' | grep MLD_MKL_LIB
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

function f_stuff_gnu_is246206 {

  # to stop stuff sometimes
  envs --create_setenv="${HOME}/mld_setenv_gnu.tmp" # get current bash environ status
  export mld_do_python_tests=ON   # to do python integration tests
  export mld_do_make=ON           # to compile milady

  export MLD_OPT_TRAC=OFF          # to log cmake all variables contents

  # check bash environment is clean and previously coherent
  [ -z ${MLD_TYPE} ] && export MLD_TYPE=GNU  # mandatory else risky
  [ ! "${MLD_TYPE}" == "GNU" ] && in_red_pause "ERROR: bash environment is not clean, continue at your risks"

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
  export MLD_OAP_ROODIR="/export/home/catA/intel"  # needed for mkl
  export MLD_OAP_INSDIR="${MLD_OAP_ROODIR}/oneapi"
  export MLD_ROODIR="${LOC_ROODIR}/MLD"
  export MLD_SRCDIR="${MLD_ROODIR}/MILADY"
  export MLD_BUIDIR="${MLD_ROODIR}/mld_build"     # is user choice
  export MLD_INSDIR="${MLD_ROODIR}/mld_install"   # is user choice
  export MLD_TESDIR="${MLD_ROODIR}/mld_tests"     # is user choice
  export MLD_MPI_ROODIR="${MLD_ROODIR}/mpi"
  export MLD_MPI_INSDIR="${MLD_MPI_ROODIR}/install_412_gnu"  # /usr/local/iopenmpi.4.1.2/lib64/
  export MLD_CMAKE_BUILD_TYPE=${MLD_CMAKE_BUILD_TYPE:-Release}  # or "Debug" #
  export MLD_SETENV_DONE=ON
  export MLD_NPROC=${MLD_NPROC:-8}

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

  if [ -z ${SETVARS_COMPLETED} ]; then
    # Intel Trace Analyzer and Collector (ITAC) is an MPI profiling and tracing tool
    in_red 'need launched setvars.sh oneAPI, set it for GNU - mkl only'
    export CONFIG_SETVARS=${MLD_ROODIR}/config_setvars_oneapi.tmp
    cat <<EOT > ${CONFIG_SETVARS}
default=exclude
mkl=latest
EOT
    in_green 'config oneAPI setvars on ${MLD_ROODIR}/config_setvars_oneapi.tmp'
    cat ${CONFIG_SETVARS} # to verify
    source ${MLD_OAP_INSDIR}/setvars.sh --config=${CONFIG_SETVARS}
  else
    in_green 'GNU compiler use, hope it is OK (without intel mpi)'
  fi

  # some verifications
  # hope there is NOT intel mpi
  envs -g '/mpi' || n_red_pause 'problem: needed there is NOT intel mpi'

  if [ -z ${MKLROOT} ]; then
    in_red_pause 'problem: inexisting ${MKLROOT} from oneAPI\nfix it.'
    return
  fi
  export MKL_ROOT=${MKLROOT}

  # set by oneAPI setvars.sh
  envs -p MKL
  envs -g oneapi/compiler \
    && in_red_pause 'KO there is intel compilers activated, tricky, fix it.' \
    || in_green 'OK there is NO intel compilers activated, expected.' \

  in_green 'set PATH & LD_LIBRARY_PATH on openmpi'
  export PATH=${MLD_MPI_INSDIR}/bin:${PATH}
  export LD_LIBRARY_PATH=${MLD_MPI_INSDIR}/lib:${LD_LIBRARY_PATH}


  tmp=$(which cc 2> /dev/null) && in_green 'cc at '${tmp} || in_red_pause 'problem no reached cc compiler'
  tmp=$(which gcc 2> /dev/null) && in_green 'gcc at '${tmp} || in_red_pause 'problem no reached gcc compiler'
  tmp=$(which gfortran 2> /dev/null) && in_green 'gfortran at '${tmp} || in_red_pause 'problem no reached gfortran compiler'
  tmp=$(which mpifort 2> /dev/null) && in_green 'mpifort at '${tmp} || in_red_pause 'problem reached useless mpifort compiler'
  tmp=$(cc --version | head -1)
  in_green "GNU cc version is ${tmp}"

  cd ${MLD_SRCDIR}/scripts
  source ./compile_milady.bash

  envs -e PATH,LD_LIBRARY_PATH

  # unset MLD_NPROC # let it automatic
  # export MLD_NPROC=10
  # export MLD_OPT_TRAC=ON
  f_compile_milady_gnu # do cmake

  if [ "${mld_do_make}" == "OFF" ]; then  # stop here ?
    cd ${MLD_ROODIR}
    return
  fi
  f_make # do make
  in_green "make install"
  make install

  in_green 'set ${MILADYPY_ROOT_DIR} on miladypy'
  export MILADYPY_ROOT_DIR=${MLD_SRCDIR}


  in_green "create_diff_setenv from "${HOME}/mld_setenv_gnu.tmp
  envs --create_diff_setenv="${HOME}/mld_setenv_gnu.tmp" # get current bash environ status


  # example find libmpi_mpifh.so
  # in_green "find libmpi_mpifh.so"
  # locate libmpi_mpifh.so
  # ldconfig -p | grep libmpi_mpifh.so

  # f_explore_milady   # long listing only for debug

  # --------------------------------------------------------------------------
  # No network interfaces were found for out-of-band communications. We require
  # at least one available network for out-of-band messaging.
  # --------------------------------------------------------------------------
  # if fedora 34 at home out of network ... ! this solve problem Ooops
  # echo 'oob_tcp_if_include = lo' >> ./install_412_gnu/etc/openmpi-mca-params.conf
  mpirun -n 1 pwd \
    && in_green 'OK for mpirun' \
    || in_red_pause 'OK for mpirun. try: "mpirun -n 1 pwd" and fix problem.'
  # try: echo "oob_tcp_if_include = lo" >> ${MLD_MPI_INSDIR}/etc/openmpi-mca-params.conf

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
  mld_verboseTest=OFF
  f_python_tests_milady
}


in_green 'to compile milady and execute tests, type:\n  f_stuff_gnu_is246206'
[ "${USER}" == "wambeke" ] && f_stuff_gnu_is246206   # do it for me
