

### modules oneAPI

To get environment and fonctionalities or libraries 
from all (...or some) products of oneAPI

- **In this previous case** See *Use a config file for setvars* (below).

- **OR** use `module load` command

**WARNING:** `module` is a command OS-linux-distribution dependent. 
May be not installed/effective on your system by default
(as Cosmin linux `suze`, where `module` command is not installed.


#### Use a config file for setvars.sh

This advantageously replaces oneAPI `module load` procedure(s).

```
export MLD_OAP_ROODIR=.../intel
export MLD_OAP_INSDIR=${MLD_OAP_ROODIR}/oneapi        # oneAPI install
export MLD_ROODIR=$(pwd)  # for example...

# case Cosmin linux suze have not module command
# source ${MLD_OAP_INSDIR}/setvars.sh on config file

export CONFIG_SETVARS=${MLD_ROODIR}/config_setvars_oneapi.tmp
cat <<EOT > ${CONFIG_SETVARS}
default=exclude
compiler=latest
mkl=latest
EOT

# avoid mpi setvars if milady_MIX compilation (without mpi intel)
# echo "mpi=latest" >> ${CONFIG_SETVARS}

cat ${CONFIG_SETVARS} # to verify

# may be shit happens, fix as restart on new/clean bash
# WARNING: setvars.sh has already been run. Skipping re-execution.
# To force a re-execution of setvars.sh, use the '--force' option.
# Using '--force' can result in excessive use of your environment variables.

source ${MLD_OAP_INSDIR}/setvars.sh --config=${CONFIG_SETVARS}

# suze have not module: this is obsolete
# module use ${MLD_OAP_INSDIR}/modulefiles
# module avail
# module load compiler/latest
# module load mkl/latest

# some verifications
envs -g mpi
envs -g open
envs -g mkl
which ifort
which icc
```


#### create all oneAPI module files from installation

Example (*use at your own risk*)

```
# https://www.rcac.purdue.edu/knowledge/brown/compile?all=true

bash

# is221713
export MLD_OAP_ROODIR=/volatile2/intel
export MLD_ROODIR=/volatile2/MILADY/MILADY_SANDBOX/wambeke

# is232972
export MLD_OAP_ROODIR=/volatile2/catA/wambeke/oneAPI
export MLD_ROODIR=/home/catA/${USER}/MLD

export MLD_OAP_INSDIR=${MLD_OAP_ROODIR}/oneapi        # oneAPI install


cd ${MLD_OAP_INSDIR}
# to create modulefiles
source modulefiles-setup.sh

# verification
ls ./modulefiles
  advisor   compiler32     debugger       dnnl-cpu-iomp  icc          intel_ipp_ia32     mkl32    tbb32
  ccl       compiler-rt    dev-utilities  dnnl-cpu-tbb   icc32        intel_ipp_intel64  mpi      vpl
  clck      compiler-rt32  dnnl           dpct           init_opencl  itac               oclfpga  vtune
  compiler  dal            dnnl-cpu-gomp  dpl            inspector    mkl                tbb

# add path
module use ${MLD_OAP_INSDIR}/modulefiles

module avail
  --------------------- /volatile2/intel/oneapi/modulefiles 
  advisor/2022.0.0           compiler-rt/latest         dnnl-cpu-gomp/2022.0.1     icc/latest                 itac/2021.5.0              tbb/latest
  advisor/latest             compiler-rt32/2022.0.1     dnnl-cpu-gomp/latest       icc32/2022.0.1             itac/latest                tbb32/2021.5.0
  ccl/2021.5.0               compiler-rt32/latest       dnnl-cpu-iomp/2022.0.1     icc32/latest               mkl/2022.0.1               tbb32/latest
  ccl/latest                 dal/2021.5.1               dnnl-cpu-iomp/latest       init_opencl/2022.0.1       mkl/latest                 vpl/2022.0.0
  clck/2021.5.0              dal/latest                 dnnl-cpu-tbb/2022.0.1      init_opencl/latest         mkl32/2022.0.1             vpl/latest
  clck/latest                debugger/2021.5.0          dnnl-cpu-tbb/latest        inspector/2022.0.0         mkl32/latest               vtune/2022.0.0
  compiler/2022.0.1          debugger/latest            dpct/2022.0.0              inspector/latest           mpi/2021.5.0               vtune/latest
  compiler/latest            dev-utilities/2021.5.1     dpct/latest                intel_ipp_ia32/2021.5.1    mpi/latest
  compiler32/2022.0.1        dev-utilities/latest       dpl/2021.6.0               intel_ipp_ia32/latest      oclfpga/2022.0.1
  compiler32/latest          dnnl/2022.0.1              dpl/latest                 intel_ipp_intel64/2021.5.1 oclfpga/latest
  compiler-rt/2022.0.1       dnnl/latest                icc/2022.0.1               intel_ipp_intel64/latest   tbb/2021.5.0

  -------------------- /usr/share/Modules/modulefiles 
  dot         module-git  module-info modules     null        use.own

  -------------------- /etc/modulefiles 
  mpi/mpich-3.0-x86_64 mpi/mpich-3.2-x86_64 mpi/mpich-x86_64     mpi/openmpi-x86_64


module whatis compiler/latest
      Load "debugger" to debug DPC++ applications with the gdb-oneapi debugger.
      Load "dpl" for additional DPC++ APIs: https://github.com/oneapi-src/oneDPL
    ------------------------------- /volatile2/catA/wambeke/oneAPI/oneapi/modulefiles --------------------------------
         compiler/latest: Configure for use with Intel 64-bit compiler(s).


module load compiler/latest
which icc
  /volatile2/catA/wambeke/oneAPI/oneapi/compiler/2022.0.1/linux/bin/intel64/icc
which ifort
  /volatile2/catA/wambeke/oneAPI/oneapi/compiler/2022.0.1/linux/bin/intel64/ifort
which mpiifort
  which: no mpiifort in (...
module load mpi/latest
  Loading mpi version 2021.5.0
which mpiifort
  /volatile2/catA/wambeke/oneAPI/oneapi/mpi/2021.5.0/bin/mpiifort

module load mkl/latest
envs -p MKL
    MKLROOT                        = /volatile2/catA/wambeke/oneAPI/oneapi/mkl/2022.0.1

etc...
```


### MIX compile intel with mkl intel and openmpi 412

**WARNING:** openmpi 412 was previously compiled with intel compilers

```
# https://cirrus.readthedocs.io/en/master/software-libraries/intel_mkl.html

export MLD_MPI_ROODIR=${MLD_ROODIR}/mpi/install_412

# in place of module load 
# TODOcvw verify useful ?
export PATH=${MLD_MPI_ROODIR}/bin:${PATH}
export LD_LIBRARY_PATH=${MLD_MPI_ROODIR}/lib:${LD_LIBRARY_PATH}
export LIBRARY_PATH=${MLD_MPI_ROODIR}/lib:${LIBRARY_PATH}

which mpifort
    /volatile2/MILADY/MILADY_SANDBOX/wambeke/mpi/install_412/bin/mpifort
which mpirun
    /volatile2/MILADY/MILADY_SANDBOX/wambeke/mpi/install_412/bin/mpirun

cd ${MLD_OAP_ROODIR}/oneapi
# source modulefiles-setup.sh  # done only one time
module use ${MLD_OAP_ROODIR}/oneapi/modulefiles
module avail

# this is load intel compilers
module load mkl/latest    # as ifort is known 

envs -p MKL
    MKLROOT=/volatile2/intel/oneapi/mkl/2022.0.1

mkdir ${MLD_ROODIR}/test
cd ${MLD_ROODIR}/test

export MLD_SRCDIR=${MLD_ROODIR}/MILADY
export TMP=${MLD_SRCDIR}/tests

export afile=test_010_hello_world
export afile=test_010_hello_world_openmpi
export afile=test_110_test_mkl_dgemm
export afile=test_115_test_mkl_zdotc
export afile=test_115_test_mkl_zdotc_openmpi

export src=${TMP}/${afile}.f90
export exe=${afile}.exe

find ${MLD_MPI_ROODIR} -name "*.mod"

# export FC=$(which gfortran)  # useless ?
man mpifort
mpifort --showme:compile
which ifort    # is not gfortran !!! 

# This will build against the serial version of MKL
rm ${exe}
# if mpi
mpifort -o ${exe} -L${MKLROOT}/lib/intel64 -lmkl_gf_lp64 -lmkl_core -lmkl_sequential ${src}
# if not mpi
ifort -o ${exe} -L${MKLROOT}/lib/intel64 -lmkl_gf_lp64 -lmkl_core -lmkl_sequential ${src}

ldd ${exe}
${exe}
mpirun -n 4 ${exe}

# works with ifort and opempi 412
mpifort -o ${exe} -L${MKLROOT}/lib/intel64 -lmkl_gf_lp64 -lmkl_core -lmkl_sequential ${src}
mpirun -n 4 ${exe}
        hello world !
        mpi_comm_world          0
        hello world ! from rank  0 of  4
        hello world !
        mpi_comm_world          0
        hello world ! from rank  1 of  4
        hello world !
        mpi_comm_world          0
        hello world !
        mpi_comm_world          0
        hello world ! from rank  3 of  4
        hello world ! from rank  2 of  4

ldd ${exe}
        linux-vdso.so.1 =>  (0x00007fff1e3b3000)
        libmkl_gf_lp64.so.2 => /volatile2/intel/oneapi/mkl/2022.0.1/lib/intel64/libmkl_gf_lp64.so.2 (0x00007f1395549000)
        libmkl_core.so.2 => /volatile2/intel/oneapi/mkl/2022.0.1/lib/intel64/libmkl_core.so.2 (0x00007f1391194000)
        libmkl_sequential.so.2 => /volatile2/intel/oneapi/mkl/2022.0.1/lib/intel64/libmkl_sequential.so.2 (0x00007f138f77a000)
        libmpi_usempif08.so.40 => /volatile2/MILADY/MILADY_SANDBOX/wambeke/mpi/install_412/lib/libmpi_usempif08.so.40 (0x00007f138f542000)
        libmpi_usempi_ignore_tkr.so.40 => /volatile2/MILADY/MILADY_SANDBOX/wambeke/mpi/install_412/lib/libmpi_usempi_ignore_tkr.so.40 (0x00007f138f336000)
        libmpi_mpifh.so.40 => /volatile2/MILADY/MILADY_SANDBOX/wambeke/mpi/install_412/lib/libmpi_mpifh.so.40 (0x00007f138f0c7000)
        libmpi.so.40 => /volatile2/MILADY/MILADY_SANDBOX/wambeke/mpi/install_412/lib/libmpi.so.40 (0x00007f138ed69000)
        libm.so.6 => /lib64/libm.so.6 (0x00007f138ea67000)
        libpthread.so.0 => /lib64/libpthread.so.0 (0x00007f138e84b000)
        libc.so.6 => /lib64/libc.so.6 (0x00007f138e47d000)
        /lib64/ld-linux-x86-64.so.2 (0x00007f13963e7000)
        libgcc_s.so.1 => /lib64/libgcc_s.so.1 (0x00007f138e267000)
        libdl.so.2 => /lib64/libdl.so.2 (0x00007f138e063000)
        libopen-rte.so.40 => /volatile2/MILADY/MILADY_SANDBOX/wambeke/mpi/install_412/lib/libopen-rte.so.40 (0x00007f138dd8f000)
        libopen-pal.so.40 => /volatile2/MILADY/MILADY_SANDBOX/wambeke/mpi/install_412/lib/libopen-pal.so.40 (0x00007f138da40000)
        librt.so.1 => /lib64/librt.so.1 (0x00007f138d838000)
        libz.so.1 => /lib64/libz.so.1 (0x00007f138d622000)
        libifport.so.5 => /volatile2/intel/oneapi/compiler/2022.0.1/linux/compiler/lib/intel64_lin/libifport.so.5 (0x00007f138d3f4000)
        libifcoremt.so.5 => /volatile2/intel/oneapi/compiler/2022.0.1/linux/compiler/lib/intel64_lin/libifcoremt.so.5 (0x00007f1396428000)
        libimf.so => /volatile2/intel/oneapi/compiler/2022.0.1/linux/compiler/lib/intel64_lin/libimf.so (0x00007f138cd66000)
        libintlc.so.5 => /volatile2/intel/oneapi/compiler/2022.0.1/linux/compiler/lib/intel64_lin/libintlc.so.5 (0x00007f138caee000)
        libsvml.so => /volatile2/intel/oneapi/compiler/2022.0.1/linux/compiler/lib/intel64_lin/libsvml.so (0x00007f138aa5f000)
        libirng.so => /volatile2/intel/oneapi/compiler/2022.0.1/linux/compiler/lib/intel64_lin/libirng.so (0x00007f138a6f5000)

! This build against the threaded version
gfortran -o ${exe} -L${MKLROOT}/lib/intel64 -lmkl_gf_lp64 -lmkl_core -lmkl_gnu_thread ${src}
ldd ${exe}
${exe}

cd ${MLD_ROODIR}/MILADY/scripts
source compile_milady.bash
etc...
```


### example set all environ for oneapi compilers

May be it is useful

```
cd  /volatile2/intel/oneapi
find . -name "*.sh" | grep compiler
  ./compiler/2022.0.1/linux/lib/oclfpga/board/intel_a10gx_pac/linux64/libexec/sign_aocx.sh
  ./compiler/2022.0.1/linux/lib/oclfpga/board/intel_a10gx_pac/linux64/libexec/setup_permissions.sh
  ./compiler/2022.0.1/linux/lib/oclfpga/board/intel_a10gx_pac/linux64/libexec/pac_bsp_init.sh
  ./compiler/2022.0.1/linux/lib/oclfpga/board/intel_s10sx_pac/linux64/libexec/sign_aocx.sh
  ./compiler/2022.0.1/linux/lib/oclfpga/board/intel_s10sx_pac/linux64/libexec/setup_permissions.sh
  ./compiler/2022.0.1/linux/lib/oclfpga/board/intel_s10sx_pac/linux64/libexec/pac_bsp_init.sh
  ./compiler/2022.0.1/linux/lib/oclfpga/fpgavars.sh
  ./compiler/2022.0.1/linux/lib/oclfpga/init_opencl.sh
  ./compiler/2022.0.1/linux/lib/oclfpga/fpga_sys_check.sh
  ./compiler/2022.0.1/sys_check/sys_check.sh
  ./compiler/2022.0.1/env/vars.sh

which ifort
  which: no ifort in (....etc)

source ./compiler/2022.0.1/env/vars.sh

which ifort
  /volatile2/intel/oneapi/compiler/2022.0.1/linux/bin/intel64/ifort
which icpc
  /volatile2/intel/oneapi/compiler/2022.0.1/linux/bin/intel64/icpc

envs -g oneapi
  CMAKE_PREFIX_PATH       contains 'oneapi' as '/volatile2/intel/oneapi/compiler/2022.0.1/linux/IntelDPCPP'
  CMPLR_ROOT              contains 'oneapi' as '/volatile2/intel/oneapi/compiler/2022.0.1'
  FPGA_VARS_DIR           contains 'oneapi' as '/volatile2/intel/oneapi/compiler/2022.0.1/linux/lib/oclfpga'
  INTELFPGAOCLSDKROOT     contains 'oneapi' as '/volatile2/intel/oneapi/compiler/2022.0.1/linux/lib/oclfpga'
  LD_LIBRARY_PATH         contains 'oneapi' as '/volatile2/intel/oneapi/compiler/2022.0.1/linux/lib'
  LD_LIBRARY_PATH         contains 'oneapi' as '/volatile2/intel/oneapi/compiler/2022.0.1/linux/lib/x64'
  LD_LIBRARY_PATH         contains 'oneapi' as '/volatile2/intel/oneapi/compiler/2022.0.1/linux/lib/oclfpga/host/linux64/lib'
  LD_LIBRARY_PATH         contains 'oneapi' as '/volatile2/intel/oneapi/compiler/2022.0.1/linux/compiler/lib/intel64_lin'
  LIBRARY_PATH            contains 'oneapi' as '/volatile2/intel/oneapi/compiler/2022.0.1/linux/compiler/lib/intel64_lin'
  LIBRARY_PATH            contains 'oneapi' as '/volatile2/intel/oneapi/compiler/2022.0.1/linux/lib'
  MANPATH                 contains 'oneapi' as '/volatile2/intel/oneapi/compiler/2022.0.1/documentation/en/man/common'
  MODULEPATH              contains 'oneapi' as '/volatile2/intel/oneapi/modulefiles'
  NLSPATH                 contains 'oneapi' as '/volatile2/intel/oneapi/compiler/2022.0.1/linux/compiler/lib/intel64_lin/locale/%l_%t/%N'
  OCL_ICD_FILENAMES       contains 'oneapi' as '/volatile2/intel/oneapi/compiler/2022.0.1/linux/lib/x64/libintelocl.so'
  PATH                    contains 'oneapi' as '/volatile2/intel/oneapi/compiler/2022.0.1/linux/lib/oclfpga/bin'
  PATH                    contains 'oneapi' as '/volatile2/intel/oneapi/compiler/2022.0.1/linux/bin/intel64'
  PATH                    contains 'oneapi' as '/volatile2/intel/oneapi/compiler/2022.0.1/linux/bin'
  PKG_CONFIG_PATH         contains 'oneapi' as '/volatile2/intel/oneapi/compiler/2022.0.1/lib/pkgconfig'

envs -g mkl
  No env_var with any 'mkl' in value

source setvars.sh 
  :: initializing oneAPI environment ...
     bash: BASH_VERSION = 4.2.46(2)-release
     args: Using "$@" for setvars.sh arguments: 
  :: advisor -- latest
  :: ccl -- latest
  :: clck -- latest
  :: compiler -- latest
  :: dal -- latest
  :: debugger -- latest
  :: dev-utilities -- latest
  :: dnnl -- latest
  :: dpcpp-ct -- latest
  :: dpl -- latest
  :: inspector -- latest
  :: ipp -- latest
  :: itac -- latest
  :: mkl -- latest
  :: mpi -- latest
  :: tbb -- latest
  :: vpl -- latest
  :: vtune -- latest
  :: oneAPI environment initialized ::
 
envs -g mkl
  CPATH                   contains 'mkl' as '/volatile2/intel/oneapi/mkl/2022.0.1/include'
  LD_LIBRARY_PATH         contains 'mkl' as '/volatile2/intel/oneapi/mkl/2022.0.1/lib/intel64'
  LIBRARY_PATH            contains 'mkl' as '/volatile2/intel/oneapi/mkl/2022.0.1/lib/intel64'
  MKLROOT                 contains 'mkl' as '/volatile2/intel/oneapi/mkl/2022.0.1'
  NLSPATH                 contains 'mkl' as '/volatile2/intel/oneapi/mkl/2022.0.1/lib/intel64/locale/%l_%t/%N'
  PATH                    contains 'mkl' as '/volatile2/intel/oneapi/mkl/2022.0.1/bin/intel64'
  PKG_CONFIG_PATH         contains 'mkl' as '/volatile2/intel/oneapi/mkl/2022.0.1/lib/pkgconfig'

```
