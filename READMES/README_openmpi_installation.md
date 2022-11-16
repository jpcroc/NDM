
## openmpi installation

Compiling with intel oneAPI compilers

### set and check intel environment (to get ifort etc)

See [README_milady_installation.md](tips/README_oneAPI_setvars_modules.md)

```
export MLD_ROODIR=.../MLD # for example
# to find envs and other scripts *.bash
export PATH=${MLD_ROODIR}/MILADY/scripts:${PATH}
envs -g open
    No env_var with any 'open' in value
cd ${MLD_ROODIR}/MILADY/scripts
ls *.bash | sort
    compile_milady.bash
    compile_openmpi.bash
    utils_milady.bash

export MLD_OAP_INSDIR=.../oneapi # for example
# see READMES/tips/README_oneAPI_setvars_modules.md
export CONFIG_SETVARS=${MLD_ROODIR}/config_setvars_oneapi.tmp
cat <<EOT > ${CONFIG_SETVARS}
default=exclude
compiler=latest
mkl=latest
EOT
cat ${CONFIG_SETVARS} # to verify

# may be...
# WARNING: setvars.sh has already been run. Skipping re-execution.
# To force a re-execution of setvars.sh, use the '--force' option.
# Using '--force' can result in excessive use of your environment variables.

source ${MLD_OAP_INSDIR}/setvars.sh --config=${CONFIG_SETVARS}

which ifort
    .../oneapi/compiler/2022.0.1/linux/bin/intel64/ifort
```

### compile openmpi with oneAPI compilers

Using [compile_openmpi.bash](MILADY/scripts/compile_openmpi.bash)

```
export MLD_ROODIR=.../MLD # for example
export MLD_SRCDIR=${MLD_ROODIR}/MILADY

# on is246206 for example
envs -p MLD_
    MLD_BUIDIR                     = /export/home/catA/wambeke/MLD/build
    MLD_INSDIR                     = /export/home/catA/wambeke/MLD/install
    MLD_ROODIR                     = /export/home/catA/wambeke/MLD
    MLD_SRCDIR                     = /export/home/catA/wambeke/MLD/MILADY
    MLD_TESDIR                     = /export/home/catA/wambeke/MLD/Tests
    
envs -e PATH
    PATH                           = /export/home/catA/wambeke/MLD/MILADY/scripts
                                     /export/home/catA/miniconda/envs/py3qt5/bin
                                     /export/home/catA/miniconda/condabin
                                     /export/home/catA/intel/oneapi/mkl/2022.0.1/bin/intel64
                                     /export/home/catA/intel/oneapi/compiler/2022.0.1/linux/lib/oclfpga/bin
                                     /export/home/catA/intel/oneapi/compiler/2022.0.1/linux/bin/intel64
                                     /export/home/catA/intel/oneapi/compiler/2022.0.1/linux/bin
                                     .
                                     /home/catA/wambeke
                                     /home/catA/wambeke/commons
                                     /home/catA/wambeke/.local/bin
                                     /home/catA/wambeke/bin
                                     /usr/lib64/qt-3.3/bin
                                     /usr/share/Modules/bin
                                     /usr/local/bin
                                     /usr/bin
                                     /bin
                                     /usr/local/sbin
                                     /usr/sbin

cd ${MLD_SRCDIR}/scripts
source ./compile_openmpi.bash

# case wambeke
f_setenv_openmpi_412_wambeke
f_wget_openmpi_412

# ...take a time...
f_compile_openmpi_intel

ls -alt ${MLD_MPI_INSDIR} # check
```