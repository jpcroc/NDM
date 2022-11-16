#!/usr/bin/env bash

# openmpi compilation etc, as bash functions
# cd .../MILADY/scripts && source ./compilation_openmpi.bash
# see https://www.open-mpi.org/faq/?category=building#easy-build

source utils_milady.bash

function f_setenv_openmpi_wambeke {
  # as van wambeke example of openmpi compilation directories: all under $MLD_MPI_ROODIR
  if [ -z "${MLD_MPI_SETENV}" ]; then  # do only one time
    # export MLD_MPI_ROODIR=/volatile2/${USER}/MLD  # MLD_MPI_ROODIR idem MLD_ROODIR by choice
    # export MLD_MPI_ROODIR=/home/catA/${USER}/MLD
    export MLD_MPI_ROODIR=${MLD_ROODIR}
    # export MLD_MPI_SRCDIR=${MLD_MPI_ROODIR}/openmpi-4.1.2
    export MLD_MPI_SRCDIR=${MLD_MPI_ROODIR}/openmpi-3.1.6
    # useless autotools export MLD_MPI_BUIDIR=${MLD_MPI_ROODIR}/mpi_build
    export MLD_MPI_INSDIR=${MLD_MPI_ROODIR}/mpi_install
    export MLD_MPI_SETENV=ON
  fi
}

function f_setenv_openmpi_412_wambeke {
  # as van wambeke example of openmpi compilation directories: all under $MLD_MPI_ROODIR
  export MLD_MPI_ROODIR=${MLD_ROODIR}/mpi
  export MLD_MPI_SRCDIR=${MLD_MPI_ROODIR}/openmpi-4.1.2
  export MLD_MPI_INSDIR=${MLD_MPI_ROODIR}/install_412
  export MLD_MPI_SETENV=ON
  envs -p MLD_MPI_
}


function f_setenv_openmpi_marinica {
  # as marinica example of openmpi compilation directories
  if [ -z "${MLD_SETENV_DONE}" ]; then  # do only one time
    in_red "f_setenv_openmpi_marinica TODO"
  fi
}

function f_setenv_openmpi {
  # for now default is as wambeke
  # could define some f_setenv for ${USER} and ${HOSTNAME} etc
  unset MLD_MPI_SETENV
  [ "${USER}" == "wambeke" ] && f_setenv_openmpi_wambeke
  [ "${USER}" == "marinica" ] && f_setenv_openmpi_marinica
  [ -z "${MLD_MPI_SETENV}" ] && in_red 'problem unknown ${USER}, may be fix your preferences in f_setenv_openmpi_'${USER}
  [ -z "${MLD_MPI_SETENV}" ] && f_setenv_openmpi_wambeke
  envs -p MLD_MPI_
}

function f_wget_openmpi_412 {
  if [ ! -d "${MLD_MPI_SRCDIR}" ]; then
    in_red "download openmpi 412 in "${MLD_MPI_ROODIR}
    [ ! -d "${MLD_MPI_ROODIR}" ] && mkdir ${MLD_MPI_ROODIR}
    cd ${MLD_MPI_ROODIR}
    wget https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.2.tar.gz
    in_red 'detar ... maybe patient...'
    tar -xf ./openmpi-4.1.2.tar.gz
    ls -alt
  fi
}

function f_wget_openmpi_316 {
  if [ ! -d "${MLD_MPI_SRCDIR}" ]; then
    in_red "download openmpi 316 in "${MLD_MPI_ROODIR}
    [ ! -d "${MLD_MPI_ROODIR}" ] && mkdir ${MLD_MPI_ROODIR}
    cd ${MLD_MPI_ROODIR}
    wget https://download.open-mpi.org/release/open-mpi/v3.1/openmpi-3.1.6.tar.gz
    in_red 'detar ...take a time, be patient...'
    tar -xf ./openmpi-3.1.6.tar.gz
    ls -alt
  fi
}

function f_explore_openmpi {
  # list libraries dependencies of executables, for verification
  in_red '\ndirectory ${MLD_MPI_SRCDIR}' ${MLD_MPI_SRCDIR}
  ls -alt ${MLD_MPI_SRCDIR}
  in_red '\ndirectory ${MLD_MPI_INSDIR}' ${MLD_MPI_INSDIR}
  tree -d ${MLD_MPI_INSDIR}
  in_red '\ndirectory ${MLD_MPI_INSDIR}/bin' ${MLD_MPI_INSDIR}/bin
  ls -alt ${MLD_MPI_INSDIR}/bin
}

function f_reset_source_openmpi {
  # get clean openmpi sources if previous configure cache pollution, or not
  if [ -z "${MLD_MPI_SRCDIR}" ]; then
    in_red_pause 'problem undefined ${MLD_MPI_SRCDIR}'
    return
  fi
  cd ${MLD_MPI_SRCDIR}
  if [ -f "Makefile" ]; then  # created by previous configure execution
    in_red_pause "there is pollution of previous compilation in ${MLD_MPI_SRCDIR}"
    # reset and recreate clean openmpi dir
    cd ${MLD_MPI_ROODIR}
    rm -rf ${MLD_MPI_SRCDIR}
    in_red_ 'supposedly openmpi tar.gz here... see f_wget_openmpi_xxx'
    [[ ${MLD_MPI_SRCDIR} = *"3.1.6"* ]] && tar -xf ./openmpi-3.1.6.tar.gz --directory ${MLD_MPI_ROODIR}
    [[ ${MLD_MPI_SRCDIR} = *"4.1.2"* ]] && tar -xf ./openmpi-4.1.2.tar.gz --directory ${MLD_MPI_ROODIR}
  else
    in_green "there is NO pollution of previous compilation in ${MLD_MPI_SRCDIR}"
  fi
}

function f_compile_openmpi {
  [ -z "${MLD_MPI_SRCDIR}" ] && in_red 'problem undefined ${MLD_MPI_SRCDIR}'
  # useless autotools [ -z "${MLD_MPI_BUIDIR}" ] && in_red 'problem undefined ${MLD_MPI_BUIDIR}'
  [ -z "${MLD_MPI_INSDIR}" ] && in_red 'problem undefined ${MLD_MPI_INSDIR}'
  [ -z "${FC}" ] && in_red 'problem undefined ${FC}'
  [ -z "${CC}" ] && in_red 'problem undefined ${CC}'
  # f_wget_openmpi_316
  MLD_NPROC=${MLD_NPROC:-8} # ${NPROC:-$(f_nproc)}

  # get clean openmpi sources if previous configure pollution
  f_reset_source_openmpi

  if [ ! -z "${MLD_MPI_SRCDIR}" ]; then
    cd ${MLD_MPI_SRCDIR}
    in_green 'FC='${FC}
    in_green 'F77='${F77}
    in_green 'CC='${CC}
    in_green 'CXX='${CXX}
    in_green 'MLD_MPI_INSDIR='${MLD_MPI_INSDIR}
    in_red 'f_compile_openmpi configure ...take a time, be patient...'
    # ./configure CC=icc CXX=icpc F77=ifort FC=ifort --prefix=${MLD_MPI_INSDIR}
    # WARNING are in stderr
    ./configure CC=${CC} CXX=${CXX} F77=${F77} FC=${FC} --prefix=${MLD_MPI_INSDIR} > ./log_openmpi_configure.tmp \
      && in_green 'OK configure see ./log_openmpi_configure.tmp' \
      || in_red_pause 'KO configure see ./log_openmpi_configure.tmp, try to detar a new/clean openmpi directory'
    in_red 'f_compile_openmpi make -j'${MLD_NPROC}' ...take a time, be patient...'
    make -j${MLD_NPROC} all > ./log_openmpi_make.tmp \
      && in_green 'OK make see ./log_openmpi_make.tmp' \
      || in_red_pause 'KO configure see ./log_openmpi_make.tmp'
    in_red 'f_compile_openmpi make install ...take a time, be patient...'
    make install > ./log_openmpi_make_install.tmp \
      && in_green 'OK make install see ./log_openmpi_make_install.tmp' \
      || in_red_pause 'KO configure see ./log_openmpi_make_install.tmp'
    in_green "End f_compile_openmpi" ${MLD_MPI_INSDIR}
  fi
}

function f_set_mpi_instdir {
  # append ${1} as "intel" "gnu" to $MLD_MPI_INSDIR
  # (storing original $MLD_MPI_INSDIR in $MLD_MPI_INSDIR_ORIG if not done)
  [ -z "${MLD_MPI_INSDIR_ORIG}" ] && export MLD_MPI_INSDIR_ORIG=${MLD_MPI_INSDIR}
  # example from ../mpi/install_412 to  .../mpi/install_412_gnu
  export MLD_MPI_INSDIR=${MLD_MPI_INSDIR_ORIG}"_"${1}
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
