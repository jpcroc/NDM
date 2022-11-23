#!/usr/bin/env bash

# ndm compilation etc, as bash functions
# cd .../NDM/scripts && source ./compile_ndm.bash

source utils_ndm.bash
# set -x for script debug

function f_setenv_ndm_wambeke {
  # as van wambeke example of ndm compilation directories: all under $NDM_ROODIR
  # export NDM_ROODIR="/volatile2/${USER}/MLD"  # something like that
  in_green "f_setenv_ndm_wambeke"
  # envs -p NDM_

  [ -z "${NDM_ROODIR}" ] && in_red_pause 'problem undefined NDM_ROODIR, define it'
  if [ -z "${NDM_SETENV_DONE}" ]; then  # do only one time
    export NDM_SRCDIR=${NDM_ROODIR}/NDM
    export NDM_BUIDIR=${NDM_ROODIR}/ndm_build
    export NDM_INSDIR=${NDM_ROODIR}/ndm_install

    # should be set before elsewhere as system configuration
    [ "${HOSTNAME}" == "is221713" ] && export NDM_OAP_ROODIR=/volatile2/intel/oneapi
    [ "${HOSTNAME}" == "is232972.intra.cea.fr" ] && export NDM_OAP_ROODIR=/volatile2/catA/wambeke/oneAPI/oneapi
    [ -z "${NDM_OAP_ROODIR}" ] && in_red 'problem undefined ${NDM_OAP_ROODIR}, define it.'
    #[ -z "${MKL_ROOT}" ] && export MKL_ROOT=${NDM_OAP_ROODIR}/mkl/latest  # latest 2022.0.1

    [ ! -d "${NDM_OAP_ROODIR}" ]  && in_red 'problem inexisting directory ${NDM_OAP_ROODIR} '${NDM_OAP_ROODIR}
    #[ ! -d "${MKL_ROOT}" ]  && in_red_pause 'problem inexisting directory ${MKL_ROOT} '${MKL_ROOT}

    # ...where you set hdf5 compilation install directory
    # (obsolete) export HDF5_ROOT=${NDM_ROODIR}/hdf_install

    # ...where you set Marinica python tests directory
    [ -z "${NDM_TESDIR}" ] && export NDM_TESDIR=${NDM_ROODIR}/ndm_tesdir
    [ "${HOSTNAME}" == "is221713" ] && export NDM_TESDIR="/volatile2/wambeke/ttmp/Tests"
    [ ! -d "${NDM_TESDIR}" ]  && in_red 'problem inexisting directory ${NDM_TESDIR} '${NDM_TESDIR}

    # ...where you set oneapi !!! only one time !!!
    if [ ! "${SETVARS_COMPLETED:-}" == "1" ]; then
      in_red 'problem oneAPI setvars not done, do it as all default mode could be dangerous for MIX'
      in_red_pause 'try: source ${NDM_OAP_ROODIR}/setvars.sh'
    fi
    export NDM_SETENV_DONE=ON
  fi
  # envs -p NDM_
}



function f_setenv_ndm {
  # for now default is as wambeke
  # could define some your customized f_setenv_xx for ${USER} and ${HOSTNAME} etc
  in_green "f_setenv_ndm"
  [ "${1}" == "FORCE" ] && unset NDM_SETENV_DONE
  [ "${USER}" == "wambeke" ] && f_setenv_ndm_wambeke
  # [ "${USER}" == "marinica" ] && f_setenv_ndm_marinica
  if [ -z "${NDM_SETENV_DONE}" ]; then
     in_red_pause 'problem unknown ${USER}, may be fix your preferences in your f_setenv_ndm_'${USER}
     f_setenv_ndm_wambeke  # as default
  fi
  envs -p NDM_  # resume
}


function f_set_ndm_buidir {
  # append ${1} as "intel" "gnu" "mix" to $NDM_BUIDIR
  # (storing original $NDM_BUIDIR in $NDM_BUIDIR_ORIG if not done)
  [ -z "${NDM_BUIDIR_ORIG}" ] && export NDM_BUIDIR_ORIG=${NDM_BUIDIR}
  # example from ../ndm_build to  .../ndm_build_mix
  [[ ! "${NDM_BUIDIR}" == *"_${1}"* ]] && export NDM_BUIDIR=${NDM_BUIDIR_ORIG}"_"${1}
  # envs -p NDM_BUIDIR
}

function f_set_ndm_insdir {
  # append ${1} as "intel" "gnu" "mix" to $NDM_INSDIR
  # (storing original $NDM_INSDIR in $NDM_INSDIR_ORIG if not done)
  [ -z "${NDM_INSDIR_ORIG}" ] && export NDM_INSDIR_ORIG=${NDM_INSDIR}
  # example from ../ndm_install to  .../ndm_install_mix
  [[ ! "${NDM_INSDIR}" == *"_${1}"* ]] && export NDM_INSDIR=${NDM_INSDIR_ORIG}"_"${1}
  # envs -p NDM_INSDIR
}

function f_clean_dir {
  # clean dir contents, no remove dir if existing, warning if $1 is file
  in_green "clean directory ${1}"
  if [ ! -z "${1}" ]; then  # avoid 'rm /*'
    [ -d "${1}" ] && rm -rf ${1}/*
    [ -f "${1}" ] && in_red "expected directory ${1} to clean, but is a file"
    [ ! -e "${1}" ] && in_green "inexisting dir ${1}, clean is superfluous"
  else
    in_red 'f_clean_dir expected directory name ($1) empty, fix it.'
  fi
}

function f_clean_ndm {
  # manual clean NDM source directory
  [ -z "${NDM_SRCDIR}" ] && in_red 'problem undefined ${NDM_SRCDIR}'
  if [ ! -z "${NDM_SRCDIR}" ]; then
    # in case of (eventually) cmake pollution in NDM source directory
    rm -rf ${NDM_SRCDIR}/mod
    rm -rf ${NDM_SRCDIR}/lib
    rm -rf ${NDM_SRCDIR}/build
    rm -rf ${NDM_SRCDIR}/CMakeFiles
    rm -f  ${NDM_SRCDIR}/CMakeCache.txt ${NDM_SRCDIR}/cmake_install.cmake
    rm -f  ${NDM_SRCDIR}/Makefile ${NDM_SRCDIR}/install_manifest.txt
    rm -f  ${NDM_SRCDIR}/libNDM.a

    # classical clean install dir all case
    [ ! -z "${NDM_INSDIR}" ] && rm -rf ${NDM_INSDIR} || in_red 'problem undefined ${NDM_INSDIR}'

    # classical clean build dir all case
    # stay bug here TODO  what where is build/mod
    if [ ! -z "${NDM_BUIDIR_ORIG}" ]; then
      [ ! "${NDM_BUIDIR_ORIG}" == "${NDM_BUIDIR}" ] && f_clean_dir ${NDM_BUIDIR_ORIG}
    fi

    if [ ! -z "${NDM_BUIDIR}" ]; then
      f_clean_dir ${NDM_BUIDIR}
    else
      in_red 'problem undefined ${NDM_BUIDIR}'
    fi
  fi
}

function f_explore_ndm {
  # list libraries dependencies of executables, for verification
  in_red '\ndirectory ${NDM_SRCDIR}' ${NDM_SRCDIR}
  ls -alt ${NDM_SRCDIR}
  in_red '\ndirectory ${NDM_BUIDIR}' ${NDM_BUIDIR}
  ls -alt ${NDM_BUIDIR}
  tree -d ${NDM_BUIDIR}
  in_red '\ndirectory ${NDM_INSDIR}' ${NDM_INSDIR}
  tree ${NDM_INSDIR}
  f_ldd_exe
}

function f_ldd_exe {
  # $1 could be -v as verbose
  tmp=$(find ${NDM_BUIDIR} -name "*.exe")
  for i in $tmp; do
    in_green '\nldd '${1} ${i}
    ldd $1 $i
  done
}

function f_end_ndm_compile {
  cd ${NDM_BUIDIR}
  in_green 'End cmake build in directory '${NDM_BUIDIR}
  ls -alt ${NDM_BUIDIR}
  in_green '\nNow you could type: TODO THIS IS NOT TESTED\n
  cmake --build . --target all --parallel '${NDM_NPROC}'\n
  cmake --build . --target install\n
  \n
  cmake --build . --target test\n
  f_ctest_ndm\n
  f_explore_ndm\n
  '
}

function f_compile_ndm {
  # do the NDM cmake compile job with remove all previous build
  in_red "f_compile_ndm TODO THIS IS NOT TESTED"

  export NDM_CMAKE_BUILD_TYPE=${NDM_CMAKE_BUILD_TYPE:-Release}  # or Debug
  # export NDM_OPT_TRAC=${NDM_OPT_TRAC:-OFF} # OFF by default as first parameter of function
  export NDM_NPROC=${NDM_NPROC:-$(f_nproc)}

  if [ -z ${NDM_TYPE} ]; then
    in_red_pause 'set ${NDM_TYPE} by default GNU'
    export NDM_TYPE=GNU
  fi

  # obsolete standart entry point for cmake find hdf5
  # [ -z "${HDF5_ROOT}" ] && in_red 'problem undefined ${HDF5_ROOT} (try ${NDM_ROODIR}/hdf_install)'

  [ -z "${NDM_SRCDIR}" ] && in_red_pause 'problem undefined ${NDM_SRCDIR}'
  [ -z "${NDM_BUIDIR}" ] && in_red_pause 'problem undefined ${NDM_BUIDIR}'
  [ -z "${NDM_INSDIR}" ] && in_red_pause 'problem undefined ${NDM_INSDIR}'
  #[ -z "${FC}" ] && in_red_pause 'problem undefined ${FC}'
  #[ -z "${CC}" ] && in_red_pause 'problem undefined ${CC}'

  [ "${NDM_TYPE}" == "GNU" ] && export NDM_PRESET_FILE=${NDM_SRCDIR}/cmake_files/ndm_preset_gnu_parallel.cmake
  [ "${NDM_TYPE}" == "INTEL" ] && export NDM_PRESET_FILE=${NDM_SRCDIR}/cmake_files/ndm_preset_oneapi_parallel.cmake

  envs -p NDM_
  # envs -g HDF_

  # for mix only export NDM_MKL_LIB='mkl_scalapack_lp64;mkl_intel_lp64;mkl_sequential;mkl_core;mkl_blacs_openmpi_lp64;pthread;iomp5;m;dl;mpi;mpi_mpifh'

  if [ ! -z "${NDM_SRCDIR}" ]; then
    f_clean_ndm
    in_green FC=${FC}
    in_green CC=${CC}
    # cmd="$(f_cmake3) -S${NDM_SRCDIR} -B${NDM_BUIDIR} -DCMAKE_BUILD_TYPE=${NDM_CMAKE_BUILD_TYPE} -DNDM_OPT_TRACE=${NDM_OPT_TRAC}"
    cmd="$(f_cmake3) -C ${NDM_PRESET_FILE} -S ${NDM_SRCDIR}"
    in_green ${cmd}
    ${cmd} && f_end_ndm_compile || in_red_pause 'problem cmake command'
    in_green "End f_compile_ndm, stay make -j${NDM_NPROC} && make install\n"
  fi
}

function f_compile_ndm_intel {
  in_green "f_compile_ndm_intel"
  in_green "compile ndm INTEL with intel mpiifort and intel mpi"
  #export FC=$(f_which mpiifort)
  #export CC=$(f_which icc)
  f_set_ndm_buidir intel
  f_set_ndm_insdir intel
  f_compile_ndm
}

function f_help_compile_gnu {
  in_red '\n
  May be you have to type:\n
  scl enable devtoolset-10 bash\n
  module load openmpi\n
  or\n
  May be gnu 9 from system:\n
  module avail gcc\n
  module load gcc/9.3.0\n
  module avail openmpi\n
  module load openmpi/gcc_9.3.0/4.0.1\n
  '
}

function f_compile_ndm_gnu {
  in_green "f_compile_ndm_gnu"
  in_green "compile ndm GNU with gnu mpifort and openmpi, this is tricky as link to oneApi mkl libraries for gnu gf"
  #export FC=$(f_which mpifort)
  #[ "${FC}" == "UNKNOWN" ] && f_help_compile_gnu
  #export CC=$(f_which cc)
  f_set_ndm_buidir gnu
  f_set_ndm_insdir gnu
  # export MKL_ROOT=/volatile2/catA/wambeke/oneAPI/oneapi/mkl/2022.0.1
  #envs -g mkl
  #[ -z ${MKLROOT} ] && in_red_pause 'problem undefined ${MK_ROOT}'
  #[ -z ${MKL_ROOT} ] && in_red_pause 'problem undefined ${MKL_ROOT}'
  #[ -z ${MKLROOT} ] && export MKLROOT=${MKL_ROOT} || export MKL_ROOT=${MKROOT}
  #[ -z ${MKL_ROOT} ] && in_red 'problem undefined ${MKL_ROOT}'
  f_compile_ndm
}


function f_make {
  [ -z "${NDM_BUIDIR}" ] && in_red 'problem undefined ${NDM_BUIDIR}'
  cd ${NDM_BUIDIR}
  in_green "make -j${NDM_NPROC}"
  make -j${NDM_NPROC}  # 10
}

function f_test_gnu_is221713 {
  # cd /volatile2/wambeke/ttmp/Tests/NDM_2020/md_cg_ml
  cd /volatile2/wambeke/ttmp/Tests/train_tests/fit08
  which mpirun
  mpirun -np 2 /volatile2/NDM/NDM_SANDBOX/wambeke/build/bin/ndm_main.exe
  cd /volatile2/NDM/NDM_SANDBOX/wambeke/build
}

function f_only_one_test {
  in_red 'only launch one test (for example)'
  in_green 'cd ${NDM_TESDIR}/NDM_2020/md_cg_ml'
  in_green 'mpirun -np 2 ${NDM_BUIDIR}//bin/ndm_main.exe'
  in_green 'cd ${NDM_BUIDIR  # as you want'
  cd ${NDM_TESDIR}/NDM_2020/md_cg_ml
  mpirun -np 2 ${NDM_BUIDIR}/bin/ndm_main.exe
}

function f_module_load_irene_intel_20 {
  in_green "*** modules load as marinica intel/20 TODO OBSOLETE***"
  module load intel/20.0.0
  module load mpi/openmpi/4.0.5
  module load mkl/20.0.0
  module load scalapack/mkl/20.0.0
  module load gnu/10.1.0
  module load python3/3.7.5 # python before cmake 3.18 as load cmake 3.15
  module load cmake/3.18.1
}

function f_module_load_irene_oneapi_21_marinica {
  in_green "*** modules load as marinica oneapi+openmpi OK 220331 ***"
  module load mpi/openmpi/4.0.5 scalapack/mkl/21.3.0
  module load fortran/inteloneapi/21.4.0
  in_green "*** modules load python 3 ***"
  module load python3/3.7.5 # python before cmake 3.18 as load cmake 3.15
  in_green "*** modules load cmake 3.18 ***"
  module load cmake/3.18.1
  in_green "*** env MKL ***"
  envs -e MKLROOT
}

function f_module_load_irene_oneapi_21_wambeke {
  in_green "*** modules load as wambeke all oneapi 220401 TODO do not works 220104 ***"
  module load inteloneapi/21.4.0
  module load mpi/intelmpi/21.4.0
  # module load mkl/21.4.0
  # in_green "*** modules load python 3 ***" conflict!
  # module load python3/3.7.5 # python before cmake 3.18 as load cmake 3.15
  in_green "*** modules load cmake 3.18 ***"
  module load cmake/3.18.1
  in_green "*** env MKL ***"
  envs -e MKLROOT
}

function f_compile_ndm_irene_intel {
  # TODO KO 220401
  export NDM_SRCDIR=${NDM_ROODIR}/NDM
  export NDM_BUIDIR=${NDM_ROODIR}/ndm_build_intel
  export NDM_INSDIR=${NDM_ROODIR}/ndm_install_intel

  # ...where you choose hdf5 compilation install directory
  # export HDF5_ROOT=${NDM_ROODIR}/hdf_install

  # ...where you put Marinica python tests directory
  export NDM_TESDIR=${NDM_ROODIR}/ndm_tesdir
  export NDM_SETENV_DONE=ON

  cd ${NDM_ROODIR}
  in_red "compile ndm irene intel with mpiifort and intelmpi"
  export FC=$(f_which mpiifort)
  [ "${FC}" == "UNKNOWN" ] && in_red "FC UNKNOWN"
  export CC=$(f_which icc)
  [ "${CC}" == "UNKNOWN" ] && in_red "CC UNKNOWN"
  # export MKL_ROOT=/volatile2/catA/wambeke/oneAPI/oneapi/mkl/2022.0.1
  #export MKLROOT=${MKL_ROOT}
  #[ -z "${MKL_ROOT}" ] && in_red 'problem undefined ${MKL_ROOT}' || f_compile_ndm
}

function f_help_ctest_ndm {
 in_red "\n
  f_help_ctest_ndm TODO 220405\n
  This help is for user local computers, NOT irene (where modules load do stuff)\n
  You have to set PATH and LD_LIBRARY_PATH before launch tests.\n
  this is tricky, may be TODO a user-defined module load file.\n
  In case of using openmpi intel_compiled (mix installation, or else gnu and intel)\n
  may be (as example)\n"

 echo -e '
  # For user local computers, NOT irene (where modules load do stuff)
  export NDM_ROODIR=/home/catA/${USER}/MLD
  # NDM_OAP_ as OneAPi
  export NDM_OAP_ROODIR=/volatile2/catA/${USER}/oneAPI
  # NDM_MPI_ as OpenMPi
  export NDM_MPI_INSDIR=${NDM_ROODIR}/mpi_install_intel
  # set LD_LIBRARY_PATH
  export LD_LIBRARY_PATH_ORIG=${LD_LIBRARY_PATH}
  export LD_LIBRARY_PATH=${NDM_OAP_ROODIR}/oneapi/mkl/2022.0.1/lib/intel64:${NDM_MPI_INSDIR}/lib:${LD_LIBRARY_PATH}
  export LD_LIBRARY_PATH=${NDM_OAP_ROODIR}/oneapi/itac/2021.5.0/bin/rtlib:${LD_LIBRARY_PATH}
  export LD_LIBRARY_PATH=${NDM_OAP_ROODIR}/oneapi/compiler/2022.0.1/linux/compiler/lib/intel64_lin:${LD_LIBRARY_PATH}
  # set mpirun PATH
  export PATH_ORIG=${PATH}
  export PATH=${NDM_MPI_INSDIR}/bin:${PATH}
  envs -e LD_LIBRARY_PATH
  envs -e PATH
  conda env list
  conda activate yourPythonEnv  # python3 with numpy, matplotlib
  f_ctest_ndm
  f_python_tests_ndm
  '
}

function f_ctest_ndm {
  # launch ctest ndm
  # could type 'f_ctest -V' as verbose
  in_green "f_ctest_ndm"
  f_help_ctest_ndm
  f_setenv_ndm
  envs -e LD_LIBRARY_PATH,PATH
  envs -g oneapi/mkl
  envs -p NDM_
  olddir=$(pwd)
  # could type $1 as f_ctest_ndm -V as verbose or f_ctest_ndm '-I 11,13'
  if [ -z "${NDM_BUIDIR}" ]; then
    in_red_pause 'problem undefined ${NDM_BUIDIR}'
  else
    cd ${NDM_BUIDIR} || in_red_pause 'problem inexisting directory ${NDM_BUIDIR} '${NDM_BUIDIR}
    ctest ${1}
  fi
  cd ${olddir}
}

function f_python_tests_ndm {
  # launch python tests ndm
  # could do 'ptes_begin=7 ptes_end=9 f_python_tests_ndm' for only test 7 to 9
  in_green "f_python_tests_ndm"
  export ndm_plotFailed=${ndm_plotFailed:-ON}    # to matplotlib display failed test s
  export ndm_verboseTest=${ndm_verboseTest:-OFF}  # to get little verbosity (avoid in production)
  python3 -c "import numpy" || in_red_pause 'Needs python3 with numpy, at least, try conda info --envs'
  [ -z "${NDM_SRCDIR}" ] && in_red_pause 'problem undefined ${NDM_SRCDIR}, fix it.\nmay be you have to call f_setenv_ndm'
  [ -z "${NDM_TESDIR}" ] && in_red_pause 'problem undefined ${NDM_TESDIR}'
  if [ ! -d "${NDM_TESDIR}" ]; then
    in_red_pause 'problem inexisting directory ${NDM_TESDIR} '${NDM_TESDIR}
  else
    olddir=$(pwd)
    cd ${NDM_TESDIR}  # /volatile2/wambeke/ttmp/Tests  # ...where you put the tests
    # in_green NDM_TESDIR=${NDM_TESDIR}
    # export NDM_TESPY=${NDM_SRCDIR}/tests/small_tests.py # obsolete
    export NDM_TESPY=${NDM_SRCDIR}/ndmpy/AllTestLauncherNDMPY.sh
    in_green 'compilation environment parameters'
    envs -p NDM_  # NDM_* are general compilation environment parameters
    in_green 'execution environment parameters'
    envs -p ndm_  # ndm_* are local execution environment parameters
    in_green "execute tests with ${NDM_TESPY}"
    in_green
    # obsolete python3 ${NDM_TESPY} clean ; python3 ${NDM_TESPY} run
    ${NDM_TESPY}  # execute python unittests with default parameters
    cd ${olddir}
  fi
}

function f_get_small_tests_ndm {
  # as small tests cosmin 220201
  [ -z "${NDM_TESDIR}" ] && in_red 'problem undefined ${NDM_TESDIR}'
  scp -r wambeke@is221713:/home/wambeke/Desktop/ndm/Tests.tgz ~/${USER}/.
  tar -xf ~/Tests.tgz
  # or tar -xf ~/Desktop/ndm/Tests.tgz
  mv Tests ${NDM_TESDIR}
}

function f_help_end_compile {
  in_green '\nNow you could type:'
  [ -z "${MODULES_COLLECTION_TARGET}" ] && MODULES_COLLECTION_TARGET="UNKNOWN"
  if [ "${MODULES_COLLECTION_TARGET}" == "irene" ]; then
    # in_green '  f_compile_ndm_irene_mix or'
    in_green '  (TODO) f_compile_ndm_irene_intel or'
    in_green '  (TODO) f_compile_ndm_irene_gnu'
  else
    in_green '  f_setenv_ndm'
    # in_green '  f_compile_ndm_mix or'
    in_green '  f_compile_ndm_intel or'
    in_green '  f_compile_ndm_gnu'
  fi

  in_green '\nAnd after build, for tests, you could type:'
  in_green '  f_ctest_ndm and f_python_tests_ndm\n'
}

source utils_ndm.bash
# f_help_end_compile
