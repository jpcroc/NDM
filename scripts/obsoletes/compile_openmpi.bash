#!/usr/bin/env bash

# openmpi compilation etc, as bash functions
# cd .../MILADY/scripts && source ./compile_openmpi.bash
# see https://www.open-mpi.org/faq/?category=building#easy-build

source utils_ndm.bash

function f_setenv_openmpi_wambeke {
  # as van wambeke example of openmpi compilation directories: all under $NDM_MPI_ROODIR
  if [ -z "${NDM_MPI_SETENV}" ]; then  # do only one time
    # export NDM_MPI_ROODIR=/volatile2/${USER}/MLD  # NDM_MPI_ROODIR idem NDM_ROODIR by choice
    # export NDM_MPI_ROODIR=/home/catA/${USER}/MLD
    export NDM_MPI_ROODIR=${NDM_ROODIR}
    # export NDM_MPI_SRCDIR=${NDM_MPI_ROODIR}/openmpi-4.1.2
    export NDM_MPI_SRCDIR=${NDM_MPI_ROODIR}/openmpi-3.1.6
    # useless autotools export NDM_MPI_BUIDIR=${NDM_MPI_ROODIR}/mpi_build
    export NDM_MPI_INSDIR=${NDM_MPI_ROODIR}/mpi_install
    export NDM_MPI_SETENV=ON
  fi
}

function f_setenv_openmpi_412_wambeke {
  # as van wambeke example of openmpi compilation directories: all under $NDM_MPI_ROODIR
  export NDM_MPI_ROODIR=${NDM_ROODIR}/mpi
  export NDM_MPI_SRCDIR=${NDM_MPI_ROODIR}/openmpi-4.1.2
  export NDM_MPI_INSDIR=${NDM_MPI_ROODIR}/install_412
  export NDM_MPI_SETENV=ON
  envs -p NDM_MPI_
}


function f_setenv_openmpi {
  # for now default is as wambeke
  # could define some f_setenv for ${USER} and ${HOSTNAME} etc
  unset NDM_MPI_SETENV
  [ "${USER}" == "wambeke" ] && f_setenv_openmpi_wambeke
  [ "${USER}" == "marinica" ] && f_setenv_openmpi_marinica
  [ -z "${NDM_MPI_SETENV}" ] && in_red 'problem unknown ${USER}, may be fix your preferences in f_setenv_openmpi_'${USER}
  [ -z "${NDM_MPI_SETENV}" ] && f_setenv_openmpi_wambeke
  envs -p NDM_MPI_
}

function f_wget_openmpi_412 {
  if [ ! -d "${NDM_MPI_SRCDIR}" ]; then
    in_red "download openmpi 412 in "${NDM_MPI_ROODIR}
    [ ! -d "${NDM_MPI_ROODIR}" ] && mkdir ${NDM_MPI_ROODIR}
    cd ${NDM_MPI_ROODIR}
    wget https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.2.tar.gz
    in_red 'detar ... maybe patient...'
    tar -xf ./openmpi-4.1.2.tar.gz
    ls -alt
  fi
}

function f_wget_openmpi_316 {
  if [ ! -d "${NDM_MPI_SRCDIR}" ]; then
    in_red "download openmpi 316 in "${NDM_MPI_ROODIR}
    [ ! -d "${NDM_MPI_ROODIR}" ] && mkdir ${NDM_MPI_ROODIR}
    cd ${NDM_MPI_ROODIR}
    wget https://download.open-mpi.org/release/open-mpi/v3.1/openmpi-3.1.6.tar.gz
    in_red 'detar ...take a time, be patient...'
    tar -xf ./openmpi-3.1.6.tar.gz
    ls -alt
  fi
}

function f_explore_openmpi {
  # list libraries dependencies of executables, for verification
  in_red '\ndirectory ${NDM_MPI_SRCDIR}' ${NDM_MPI_SRCDIR}
  ls -alt ${NDM_MPI_SRCDIR}
  in_red '\ndirectory ${NDM_MPI_INSDIR}' ${NDM_MPI_INSDIR}
  tree -d ${NDM_MPI_INSDIR}
  in_red '\ndirectory ${NDM_MPI_INSDIR}/bin' ${NDM_MPI_INSDIR}/bin
  ls -alt ${NDM_MPI_INSDIR}/bin
}

function f_reset_source_openmpi {
  # get clean openmpi sources if previous configure cache pollution, or not
  if [ -z "${NDM_MPI_SRCDIR}" ]; then
    in_red_pause 'problem undefined ${NDM_MPI_SRCDIR}'
    return
  fi
  cd ${NDM_MPI_SRCDIR}
  if [ -f "Makefile" ]; then  # created by previous configure execution
    in_red_pause "there is pollution of previous compilation in ${NDM_MPI_SRCDIR}"
    # reset and recreate clean openmpi dir
    cd ${NDM_MPI_ROODIR}
    rm -rf ${NDM_MPI_SRCDIR}
    in_red_ 'supposedly openmpi tar.gz here... see f_wget_openmpi_xxx'
    [[ ${NDM_MPI_SRCDIR} = *"3.1.6"* ]] && tar -xf ./openmpi-3.1.6.tar.gz --directory ${NDM_MPI_ROODIR}
    [[ ${NDM_MPI_SRCDIR} = *"4.1.2"* ]] && tar -xf ./openmpi-4.1.2.tar.gz --directory ${NDM_MPI_ROODIR}
  else
    in_green "there is NO pollution of previous compilation in ${NDM_MPI_SRCDIR}"
  fi
}

function f_compile_openmpi {
  [ -z "${NDM_MPI_SRCDIR}" ] && in_red 'problem undefined ${NDM_MPI_SRCDIR}'
  # useless autotools [ -z "${NDM_MPI_BUIDIR}" ] && in_red 'problem undefined ${NDM_MPI_BUIDIR}'
  [ -z "${NDM_MPI_INSDIR}" ] && in_red 'problem undefined ${NDM_MPI_INSDIR}'
  [ -z "${FC}" ] && in_red 'problem undefined ${FC}'
  [ -z "${CC}" ] && in_red 'problem undefined ${CC}'
  # f_wget_openmpi_316
  NDM_NPROC=${NDM_NPROC:-8} # ${NPROC:-$(f_nproc)}

  # get clean openmpi sources if previous configure pollution
  f_reset_source_openmpi

  if [ ! -z "${NDM_MPI_SRCDIR}" ]; then
    cd ${NDM_MPI_SRCDIR}
    in_green 'FC='${FC}
    in_green 'F77='${F77}
    in_green 'CC='${CC}
    in_green 'CXX='${CXX}
    in_green 'NDM_MPI_INSDIR='${NDM_MPI_INSDIR}
    in_red 'f_compile_openmpi configure ...take a time, be patient...'
    # ./configure CC=icc CXX=icpc F77=ifort FC=ifort --prefix=${NDM_MPI_INSDIR}
    # WARNING are in stderr
    ./configure CC=${CC} CXX=${CXX} F77=${F77} FC=${FC} --prefix=${NDM_MPI_INSDIR} > ./log_openmpi_configure.tmp \
      && in_green 'OK configure see ./log_openmpi_configure.tmp' \
      || in_red_pause 'KO configure see ./log_openmpi_configure.tmp, try to detar a new/clean openmpi directory'
    in_red 'f_compile_openmpi make -j'${NDM_NPROC}' ...take a time, be patient...'
    make -j${NDM_NPROC} all > ./log_openmpi_make.tmp \
      && in_green 'OK make see ./log_openmpi_make.tmp' \
      || in_red_pause 'KO configure see ./log_openmpi_make.tmp'
    in_red 'f_compile_openmpi make install ...take a time, be patient...'
    make install > ./log_openmpi_make_install.tmp \
      && in_green 'OK make install see ./log_openmpi_make_install.tmp' \
      || in_red_pause 'KO configure see ./log_openmpi_make_install.tmp'
    in_green "End f_compile_openmpi" ${NDM_MPI_INSDIR}
  fi
}

function f_set_mpi_instdir {
  # append ${1} as "intel" "gnu" to $NDM_MPI_INSDIR
  # (storing original $NDM_MPI_INSDIR in $NDM_MPI_INSDIR_ORIG if not done)
  [ -z "${NDM_MPI_INSDIR_ORIG}" ] && export NDM_MPI_INSDIR_ORIG=${NDM_MPI_INSDIR}
  # example from ../mpi/install_412 to  .../mpi/install_412_gnu
  export NDM_MPI_INSDIR=${NDM_MPI_INSDIR_ORIG}"_"${1}
}

function f_compile_openmpi_intel {
  in_green "compile openmpi with intel ifort"
  export FC=$(f_which ifort)
  export F77=$(f_which ifort)
  export CC=$(f_which icc)
  export CXX=$(f_which icpc)
  f_set_mpi_instdir intel
  f_compile_openmpi
}

function f_compile_openmpi_gnu {
  in_green "compile openmpi with gnu gfortran, verify environment without oneAPI"
  export FC=$(f_which gfortran)
  export F77=$(f_which gfortran)
  export CC=$(f_which cc)
  export CXX=$(f_which g++)
  f_set_mpi_instdir gnu
  envs -g oneAPI
  in_red "verify environment without oneAPI please"
  f_compile_openmpi
}

in_green "Now you could type:\n
  f_setenv_openmpi\n
  f_compile_openmpi_intel\n
  f_compile_openmpi_gnu\n
"
