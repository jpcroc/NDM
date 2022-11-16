#!/usr/bin/env bash

# milady compilation etc, as bash functions
# cd .../MILADY/scripts && source ./compilation_milady.bash

source utils_milady.bash
# set -x for script debug

function f_setenv_milady_wambeke {
  # as van wambeke example of milady compilation directories: all under $MLD_ROODIR
  # export MLD_ROODIR="/volatile2/${USER}/MLD"  # something like that
  in_green "f_setenv_milady_wambeke"
  # envs -p MLD_

  [ -z "${MLD_ROODIR}" ] && in_red_pause 'problem undefined MLD_ROODIR, define it'
  if [ -z "${MLD_SETENV_DONE}" ]; then  # do only one time
    export MLD_SRCDIR=${MLD_ROODIR}/MILADY
    export MLD_BUIDIR=${MLD_ROODIR}/mld_build
    export MLD_INSDIR=${MLD_ROODIR}/mld_install

    # should be set before elsewhere as system configuration
    [ "${HOSTNAME}" == "is221713" ] && export MLD_OAP_ROODIR=/volatile2/intel/oneapi
    [ "${HOSTNAME}" == "is232972.intra.cea.fr" ] && export MLD_OAP_ROODIR=/volatile2/catA/wambeke/oneAPI/oneapi
    [ -z "${MLD_OAP_ROODIR}" ] && in_red 'problem undefined ${MLD_OAP_ROODIR}, define it.'
    [ -z "${MKL_ROOT}" ] && export MKL_ROOT=${MLD_OAP_ROODIR}/mkl/latest  # latest 2022.0.1

    [ ! -d "${MLD_OAP_ROODIR}" ]  && in_red 'problem inexisting directory ${MLD_OAP_ROODIR} '${MLD_OAP_ROODIR}
    [ ! -d "${MKL_ROOT}" ]  && in_red_pause 'problem inexisting directory ${MKL_ROOT} '${MKL_ROOT}

    # ...where you set hdf5 compilation install directory
    # (obsolete) export HDF5_ROOT=${MLD_ROODIR}/hdf_install

    # ...where you set Marinica python tests directory
    [ -z "${MLD_TESDIR}" ] && export MLD_TESDIR=${MLD_ROODIR}/mld_tesdir
    [ "${HOSTNAME}" == "is221713" ] && export MLD_TESDIR="/volatile2/wambeke/ttmp/Tests"
    [ ! -d "${MLD_TESDIR}" ]  && in_red 'problem inexisting directory ${MLD_TESDIR} '${MLD_TESDIR}

    # ...where you set oneapi !!! only one time !!!
    if [ ! "${SETVARS_COMPLETED:-}" == "1" ]; then
      in_red 'problem oneAPI setvars not done, do it as all default mode could be dangerous for MIX'
      in_red_pause 'try: source ${MLD_OAP_ROODIR}/setvars.sh'
    fi
    export MLD_SETENV_DONE=ON
  fi
  # envs -p MLD_
}



function f_setenv_milady {
  # for now default is as wambeke
  # could define some your customized f_setenv_xx for ${USER} and ${HOSTNAME} etc
  in_green "f_setenv_milady"
  [ "${1}" == "FORCE" ] && unset MLD_SETENV_DONE
  [ "${USER}" == "wambeke" ] && f_setenv_milady_wambeke
  [ "${USER}" == "marinica" ] && f_setenv_milady_marinica
  if [ -z "${MLD_SETENV_DONE}" ]; then
     in_red_pause 'problem unknown ${USER}, may be fix your preferences in your f_setenv_milady_'${USER}
     f_setenv_milady_wambeke  # as default
  fi
  envs -p MLD_  # resume
}


function f_set_mld_buidir {
  # append ${1} as "intel" "gnu" "mix" to $MLD_BUIDIR
  # (storing original $MLD_BUIDIR in $MLD_BUIDIR_ORIG if not done)
  [ -z "${MLD_BUIDIR_ORIG}" ] && export MLD_BUIDIR_ORIG=${MLD_BUIDIR}
  # example from ../mld_build to  .../mld_build_mix
  [[ ! "${MLD_BUIDIR}" == *"_${1}"* ]] && export MLD_BUIDIR=${MLD_BUIDIR_ORIG}"_"${1}
  # envs -p MLD_BUIDIR
}

function f_set_mld_insdir {
  # append ${1} as "intel" "gnu" "mix" to $MLD_INSDIR
  # (storing original $MLD_INSDIR in $MLD_INSDIR_ORIG if not done)
  [ -z "${MLD_INSDIR_ORIG}" ] && export MLD_INSDIR_ORIG=${MLD_INSDIR}
  # example from ../mld_install to  .../mld_install_mix
  [[ ! "${MLD_INSDIR}" == *"_${1}"* ]] && export MLD_INSDIR=${MLD_INSDIR_ORIG}"_"${1}
  # envs -p MLD_INSDIR
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

function f_clean_milady {
  # manual clean MILADY source directory
  [ -z "${MLD_SRCDIR}" ] && in_red 'problem undefined ${MLD_SRCDIR}'
  if [ ! -z "${MLD_SRCDIR}" ]; then
    # in case of (eventually) cmake pollution in MILADY source directory
    rm -rf ${MLD_SRCDIR}/mod
    rm -rf ${MLD_SRCDIR}/lib
    rm -rf ${MLD_SRCDIR}/build
    rm -rf ${MLD_SRCDIR}/CMakeFiles
    rm -f  ${MLD_SRCDIR}/CMakeCache.txt ${MLD_SRCDIR}/cmake_install.cmake
    rm -f  ${MLD_SRCDIR}/Makefile ${MLD_SRCDIR}/install_manifest.txt
    rm -f  ${MLD_SRCDIR}/libMILADY.a

    # classical clean install dir all case
    [ ! -z "${MLD_INSDIR}" ] && rm -rf ${MLD_INSDIR} || in_red 'problem undefined ${MLD_INSDIR}'

    # classical clean build dir all case
    # stay bug here TODO  what where is build/mod
    if [ ! -z "${MLD_BUIDIR_ORIG}" ]; then
      [ ! "${MLD_BUIDIR_ORIG}" == "${MLD_BUIDIR}" ] && f_clean_dir ${MLD_BUIDIR_ORIG}
    fi

    if [ ! -z "${MLD_BUIDIR}" ]; then
      f_clean_dir ${MLD_BUIDIR}
    else
      in_red 'problem undefined ${MLD_BUIDIR}'
    fi
  fi
}

function f_explore_milady {
  # list libraries dependencies of executables, for verification
  in_red '\ndirectory ${MLD_SRCDIR}' ${MLD_SRCDIR}
  ls -alt ${MLD_SRCDIR}
  in_red '\ndirectory ${MLD_BUIDIR}' ${MLD_BUIDIR}
  ls -alt ${MLD_BUIDIR}
  tree -d ${MLD_BUIDIR}
  in_red '\ndirectory ${MLD_INSDIR}' ${MLD_INSDIR}
  tree ${MLD_INSDIR}
  f_ldd_exe
}

function f_ldd_exe {
  # $1 could be -v as verbose
  tmp=$(find ${MLD_BUIDIR} -name "*.exe")
  for i in $tmp; do
    in_green '\nldd '${1} ${i}
    ldd $1 $i
  done
}

function f_end_mld_compile {
  cd ${MLD_BUIDIR}
  in_green 'End cmake build in directory '${MLD_BUIDIR}
  ls -alt ${MLD_BUIDIR}
  in_green '\nNow you could type:\n
  make -j'${MLD_NPROC}'\n
  make install\n
  f_ctest_milady\n
  f_explore_milady\n
  '
}

function f_compile_milady {
  # do the MILADY cmake compile job with remove all previous build
  in_green "f_compile_milady"

  export MLD_CMAKE_BUILD_TYPE=${MLD_CMAKE_BUILD_TYPE:-Release}  # or Debug
  # export MLD_OPT_TRAC=${MLD_OPT_TRAC:-OFF} # OFF by default as first parameter of function
  export MLD_NPROC=${MLD_NPROC:-$(f_nproc)}

  if [ -z ${MLD_TYPE} ]; then
    in_red_pause 'set ${MLD_TYPE} by default MIX'
    export MLD_TYPE=MIX
  fi

  # obsolete standart entry point for cmake find hdf5
  # [ -z "${HDF5_ROOT}" ] && in_red 'problem undefined ${HDF5_ROOT} (try ${MLD_ROODIR}/hdf_install)'

  [ -z "${MLD_SRCDIR}" ] && in_red_pause 'problem undefined ${MLD_SRCDIR}'
  [ -z "${MLD_BUIDIR}" ] && in_red_pause 'problem undefined ${MLD_BUIDIR}'
  [ -z "${MLD_INSDIR}" ] && in_red_pause 'problem undefined ${MLD_INSDIR}'
  [ -z "${FC}" ] && in_red_pause 'problem undefined ${FC}'
  [ -z "${CC}" ] && in_red_pause 'problem undefined ${CC}'

  envs -p MLD_
  # envs -g HDF_

  # for mix only export MLD_MKL_LIB='mkl_scalapack_lp64;mkl_intel_lp64;mkl_sequential;mkl_core;mkl_blacs_openmpi_lp64;pthread;iomp5;m;dl;mpi;mpi_mpifh'

  if [ ! -z "${MLD_SRCDIR}" ]; then
    f_clean_milady
    in_green FC=${FC}
    in_green CC=${CC}
    cmd="$(f_cmake3) -S${MLD_SRCDIR} -B${MLD_BUIDIR} -DCMAKE_BUILD_TYPE=${MLD_CMAKE_BUILD_TYPE} -DMLD_OPT_TRACE=${MLD_OPT_TRAC}"
    in_green ${cmd}
    ${cmd} && f_end_mld_compile || in_red_pause 'problem cmake command'
    in_green "End f_compile_milady, stay make -j${MLD_NPROC} && make install\n"
  fi
}

function f_compile_milady_intel {
  in_green "f_compile_milady_intel"
  in_green "compile milady INTEL with intel mpiifort and intel mpi"
  export FC=$(f_which mpiifort)
  export CC=$(f_which icc)
  f_set_mld_buidir intel
  f_set_mld_insdir intel
  f_compile_milady
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

function f_compile_milady_gnu {
  in_green "f_compile_milady_gnu"
  in_green "compile milady GNU with gnu mpifort and openmpi, this is tricky as link to oneApi mkl libraries for gnu gf"
  export FC=$(f_which mpifort)
  [ "${FC}" == "UNKNOWN" ] && f_help_compile_gnu
  export CC=$(f_which cc)
  f_set_mld_buidir gnu
  f_set_mld_insdir gnu
  # export MKL_ROOT=/volatile2/catA/wambeke/oneAPI/oneapi/mkl/2022.0.1
  envs -g mkl
  [ -z ${MKLROOT} ] && in_red_pause 'problem undefined ${MK_ROOT}'
  [ -z ${MKL_ROOT} ] && in_red_pause 'problem undefined ${MKL_ROOT}'
  #[ -z ${MKLROOT} ] && export MKLROOT=${MKL_ROOT} || export MKL_ROOT=${MKROOT}
  #[ -z ${MKL_ROOT} ] && in_red 'problem undefined ${MKL_ROOT}'
  f_compile_milady
}

function f_compile_milady_mix {
  in_green "f_compile_milady_mix"
  in_green "compile milady MIX with intel ifort and intel-compiled openmpi"
  # theses env vars could be used for implicit cmake find package
  [ -z "${MKLROOT}" ] && in_red_pause 'problem undefined ${MKLROOT}'
  export MKL_ROOT=${MKLROOT}
  [ -z "${MLD_MPI_INSDIR}" ] && in_red_pause 'problem undefined ${MLD_MPI_INSDIR}'
  export FC=$(f_which ifort)
  export CC=$(f_which icc)
  f_set_mld_buidir mix
  f_set_mld_insdir mix
  export MPI_HOME=${MLD_MPI_INSDIR}
  export MLD_MPI_ROOT=${MLD_MPI_INSDIR}    # /usr/local/iopenmpi.4.1.2/lib64/
  # in_red 'TODO have to set other way ${MLD_MPI_ROOT} ? ' ${MLD_MPI_ROOT}
  f_compile_milady
}

function f_make {
  [ -z "${MLD_BUIDIR}" ] && in_red 'problem undefined ${MLD_BUIDIR}'
  cd ${MLD_BUIDIR}
  in_green "make -j${MLD_NPROC}"
  make -j${MLD_NPROC}  # 10
}

function f_test_gnu_is221713 {
  # cd /volatile2/wambeke/ttmp/Tests/NDM_2020/md_cg_ml
  cd /volatile2/wambeke/ttmp/Tests/train_tests/fit08
  which mpirun
  mpirun -np 2 /volatile2/MILADY/MILADY_SANDBOX/wambeke/build/bin/milady_main.exe
  cd /volatile2/MILADY/MILADY_SANDBOX/wambeke/build
}

function f_only_one_test {
  in_red 'only launch one test (for example)'
  in_green 'cd ${MLD_TESDIR}/NDM_2020/md_cg_ml'
  in_green 'mpirun -np 2 ${MLD_BUIDIR}//bin/milady_main.exe'
  in_green 'cd ${MLD_BUIDIR  # as you want'
  cd ${MLD_TESDIR}/NDM_2020/md_cg_ml
  mpirun -np 2 ${MLD_BUIDIR}/bin/milady_main.exe
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

function f_compile_milady_irene_mix {
  # OK 220401
  export MLD_SRCDIR=${MLD_ROODIR}/MILADY
  export MLD_BUIDIR=${MLD_ROODIR}/mld_build_mix
  export MLD_INSDIR=${MLD_ROODIR}/mld_install_mix

  # ...where you choose hdf5 compilation install directory
  # export HDF5_ROOT=${MLD_ROODIR}/hdf_install

  # ...where you put Marinica python tests directory
  export MLD_TESDIR=${MLD_ROODIR}/mld_tesdir
  export MLD_SETENV_DONE=ON

  cd ${MLD_ROODIR}
  in_red "compile milady irene mix with ifort and openmpi"
  export FC=$(f_which ifort)
  [ "${FC}" == "UNKNOWN" ] && in_red "FC UNKNOWN"
  export CC=$(f_which icc)
  [ "${CC}" == "UNKNOWN" ] && in_red "CC UNKNOWN"

  # export MKL_ROOT=/volatile2/catA/wambeke/oneAPI/oneapi/mkl/2022.0.1
  export MKLROOT=${MKL_ROOT}
  [ -z "${MKL_ROOT}" ] && in_red 'problem undefined ${MKL_ROOT}' || f_compile_milady
}

function f_compile_milady_irene_intel {
  # TODO KO 220401
  export MLD_SRCDIR=${MLD_ROODIR}/MILADY
  export MLD_BUIDIR=${MLD_ROODIR}/mld_build_intel
  export MLD_INSDIR=${MLD_ROODIR}/mld_install_intel

  # ...where you choose hdf5 compilation install directory
  # export HDF5_ROOT=${MLD_ROODIR}/hdf_install

  # ...where you put Marinica python tests directory
  export MLD_TESDIR=${MLD_ROODIR}/mld_tesdir
  export MLD_SETENV_DONE=ON

  cd ${MLD_ROODIR}
  in_red "compile milady irene intel with mpiifort and intelmpi"
  export FC=$(f_which mpiifort)
  [ "${FC}" == "UNKNOWN" ] && in_red "FC UNKNOWN"
  export CC=$(f_which icc)
  [ "${CC}" == "UNKNOWN" ] && in_red "CC UNKNOWN"
  # export MKL_ROOT=/volatile2/catA/wambeke/oneAPI/oneapi/mkl/2022.0.1
  export MKLROOT=${MKL_ROOT}
  [ -z "${MKL_ROOT}" ] && in_red 'problem undefined ${MKL_ROOT}' || f_compile_milady
}

function f_help_ctest_milady {
 in_red "\n
  f_help_ctest_milady TODO 220405\n
  This help is for user local computers, NOT irene (where modules load do stuff)\n
  You have to set PATH and LD_LIBRARY_PATH before launch tests.\n
  this is tricky, may be TODO a user-defined module load file.\n
  In case of using openmpi intel_compiled (mix installation, or else gnu and intel)\n
  may be (as example)\n"

 echo -e '
  # For user local computers, NOT irene (where modules load do stuff)
  export MLD_ROODIR=/home/catA/${USER}/MLD
  # MLD_OAP_ as OneAPi
  export MLD_OAP_ROODIR=/volatile2/catA/${USER}/oneAPI
  # MLD_MPI_ as OpenMPi
  export MLD_MPI_INSDIR=${MLD_ROODIR}/mpi_install_intel
  # set LD_LIBRARY_PATH
  export LD_LIBRARY_PATH_ORIG=${LD_LIBRARY_PATH}
  export LD_LIBRARY_PATH=${MLD_OAP_ROODIR}/oneapi/mkl/2022.0.1/lib/intel64:${MLD_MPI_INSDIR}/lib:${LD_LIBRARY_PATH}
  export LD_LIBRARY_PATH=${MLD_OAP_ROODIR}/oneapi/itac/2021.5.0/bin/rtlib:${LD_LIBRARY_PATH}
  export LD_LIBRARY_PATH=${MLD_OAP_ROODIR}/oneapi/compiler/2022.0.1/linux/compiler/lib/intel64_lin:${LD_LIBRARY_PATH}
  # set mpirun PATH
  export PATH_ORIG=${PATH}
  export PATH=${MLD_MPI_INSDIR}/bin:${PATH}
  envs -e LD_LIBRARY_PATH
  envs -e PATH
  conda env list
  conda activate yourPythonEnv  # python3 with numpy, matplotlib
  f_ctest_milady
  f_python_tests_milady
  '
}

function f_ctest_milady {
  # launch ctest milady
  # could type 'f_ctest -V' as verbose
  in_green "f_ctest_milady"
  f_help_ctest_milady
  f_setenv_milady
  envs -e LD_LIBRARY_PATH,PATH
  envs -g oneapi/mkl
  envs -p MLD_
  olddir=$(pwd)
  # could type $1 as f_ctest_milady -V as verbose or f_ctest_milady '-I 11,13'
  if [ -z "${MLD_BUIDIR}" ]; then
    in_red_pause 'problem undefined ${MLD_BUIDIR}'
  else
    cd ${MLD_BUIDIR} || in_red_pause 'problem inexisting directory ${MLD_BUIDIR} '${MLD_BUIDIR}
    ctest ${1}
  fi
  cd ${olddir}
}

function f_activate_python3_marinica {
  # cosmin anaconda3 installation to get python3 with numpy etc
  export MARINICA_ANACONDA="/ccc/work/cont002/den/marinica/anaconda3"
  export PATH="${MARINICA_ANACONDA}/bin:${PATH}"
  source ${WAM_HOME}/conda_activate.bash
  in_green "*** python anaconda (marinica) ***"
  which python
  # conda info --envs
  # list installed packages (as examples)
  # ls ${MARINICA_ANACONDA}/lib/python3.7/site-packages | grep Qt
  # conda list | grep numpy
  python -c "import numpy" || in_red 'problem needs python3 with numpy, at least, try conda info --envs'
}

function f_python_tests_milady {
  # launch python tests milady
  # could do 'ptes_begin=7 ptes_end=9 f_python_tests_milady' for only test 7 to 9
  in_green "f_python_tests_milady"
  export mld_plotFailed=${mld_plotFailed:-ON}    # to matplotlib display failed test s
  export mld_verboseTest=${mld_verboseTest:-OFF}  # to get little verbosity (avoid in production)
  python3 -c "import numpy" || in_red_pause 'Needs python3 with numpy, at least, try conda info --envs'
  [ -z "${MLD_SRCDIR}" ] && in_red_pause 'problem undefined ${MLD_SRCDIR}, fix it.\nmay be you have to call f_setenv_milady'
  [ -z "${MLD_TESDIR}" ] && in_red_pause 'problem undefined ${MLD_TESDIR}'
  if [ ! -d "${MLD_TESDIR}" ]; then
    in_red_pause 'problem inexisting directory ${MLD_TESDIR} '${MLD_TESDIR}
  else
    olddir=$(pwd)
    cd ${MLD_TESDIR}  # /volatile2/wambeke/ttmp/Tests  # ...where you put the tests
    # in_green MLD_TESDIR=${MLD_TESDIR}
    # export MLD_TESPY=${MLD_SRCDIR}/tests/small_tests.py # obsolete
    export MLD_TESPY=${MLD_SRCDIR}/miladypy/AllTestLauncherMILADYPY.sh
    in_green 'compilation environment parameters'
    envs -p MLD_  # MLD_* are general compilation environment parameters
    in_green 'execution environment parameters'
    envs -p mld_  # mld_* are local execution environment parameters
    in_green "execute tests with ${MLD_TESPY}"
    in_green
    # obsolete python3 ${MLD_TESPY} clean ; python3 ${MLD_TESPY} run
    ${MLD_TESPY}  # execute python unittests with default parameters
    cd ${olddir}
  fi
}

function f_get_small_tests_milady {
  # as small tests cosmin 220201
  [ -z "${MLD_TESDIR}" ] && in_red 'problem undefined ${MLD_TESDIR}'
  scp -r wambeke@is221713:/home/wambeke/Desktop/milady/Tests.tgz ~/${USER}/.
  tar -xf ~/Tests.tgz
  # or tar -xf ~/Desktop/milady/Tests.tgz
  mv Tests ${MLD_TESDIR}
}

function f_help_end_compile {
  in_green '\nNow you could type:'
  [ -z "${MODULES_COLLECTION_TARGET}" ] && MODULES_COLLECTION_TARGET="UNKNOWN"
  if [ "${MODULES_COLLECTION_TARGET}" == "irene" ]; then
    in_green '  f_compile_milady_irene_mix or'
    in_green '  (TODO) f_compile_milady_irene_intel or'
    in_green '  (TODO) f_compile_milady_irene_gnu'
  else
    in_green '  f_setenv_milady'
    in_green '  f_compile_milady_mix or'
    in_green '  f_compile_milady_intel or'
    in_green '  f_compile_milady_gnu'
  fi

  in_green '\nAnd after build, for tests, you could type:'
  in_green '  f_ctest_milady and f_python_tests_milady\n'
}

source utils_milady.bash
# f_help_end_compile
