## oneAPI installation

INSTALL oneAPI on linux from scratch 
- Without sudo root
- on local disk directory


Before oneAPI installation
- See https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html?operatingsystem=linux&distributions=webdownload&options=online

After oneAPI installation
- See https://www.intel.com/content/www/us/en/develop/documentation/get-started-with-intel-oneapi-base-linux/top.html


```
export MLD_OAP_TMPDIR=/volatile2/oneAPI  # example is221713 wambeke centos7
export MLD_OAP_TMPDIR=/volatile2/catA/wambeke/oneAPI  # example is232972.intra.cea.fr wambeke

# MLD_OAP_TMPDIR stores l_BaseKit_p_2022.1.1.119.sh  l_HPCKit_p_2022.1.1.97.sh
export MLD_OAP_TMPDIR=/export/home/catA/oneAPI  # example is246206.intra.cea.fr wambeke
export MLD_OAP_ROODIR=/export/home/catA/intel   # example is246206.intra.cea.fr wambeke

mkdir ${MLD_OAP_TMPDIR}
cd ${MLD_OAP_TMPDIR}

# c and c++
wget https://registrationcenter-download.intel.com/akdlm/irc_nas/18445/l_BaseKit_p_2022.1.1.119.sh
# and also hpc contains fortran !
wget https://registrationcenter-download.intel.com/akdlm/irc_nas/18438/l_HPCKit_p_2022.1.1.97.sh

# avoid that! (directly as root) ->  sudo l_BaseKit_p_2022.1.1.119.sh
onefile=l_BaseKit_p_2022.1.1.119.sh
twofile=l_HPCKit_p_2022.1.1.97.sh

chmod a+x ${onefile}
chmod a+x ${twofile}

function get_head {
  # extract head asci file from begin to "__CONTENT__"
  # example get_head l_BaseKit_p_2022.1.1.119.sh
  # set -x
  ends="__CONTENT__"
  lines=$1
  cat $lines | while read -r line; do
    echo "$line"
    if [ "$line" == "$ends" ]; then
      break
    fi
  done
  # set +x 
}

# show to explore... optionally
# get_head ${onefile} > ${onefile}.tmp
# less ${onefile}.tmp  # show to explore as edition

# remove previous oneAPI installations pollutions
rm -rf ~/intel
rm -rf ${MLD_OAP_ROODIR}/*

# !!!!! WARNING !!!!
# do not forget to modify install-dir (as NOT effective from CLI in GUI)
# choose custom what you want
# ALSO Cosmin tell me avoid binding oneAPI python

echo 'DO NOT FORGET to change install-dir in GUI-customize : '${MLD_OAP_ROODIR}

./${onefile} -s -a \
 --install-dir=${MLD_OAP_ROODIR} \
 --download-cache=${MLD_OAP_ROODIR}/download-cache \
 --intel-sw-improvement-program-consent=decline \
 --log-dir=${MLD_OAP_ROODIR}/oneapi_logs
 
# verify results
du -sh ~/intel ${MLD_OAP_ROODIR}/*
    123M	~/intel
    2,8G	.../intel/download-cache
    16G	    .../intel/oneapi
    25M	    .../intel/oneapi_logs

find ${MLD_OAP_ROODIR}/oneapi -name 'lib' | grep mkl
    .../oneapi/mkl/2022.0.1/lib

find ${MLD_OAP_ROODIR}/oneapi -name 'compiler'
    .../oneapi/compiler
    .../oneapi/compiler/2022.0.1/modulefiles/compiler
    .../oneapi/compiler/2022.0.1/linux/compiler
# etc... 
# /export/home/catA/intel/oneapi/compiler/2022.0.1/linux/compiler # for is246206

# now get fortran oneAPI compiler
# https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html?operatingsystem=linux&distributions=webdownload&options=online
# and also hpc contains fortran!
./${twofile}

# init environment to use oneAPI
source ${MLD_OAP_ROODIR}/oneapi/setvars.sh

> mpi[... and completion]
    mpicc          mpiexec.hydra  mpifc          mpiicc         mpirun         
    mpicxx         mpif77         mpigcc         mpiicpc        mpitune        
    mpiexec        mpif90         mpigxx         mpiifort       mpitune_fast   

which mpiifort
    .../oneapi/mpi/2021.5.0/bin/mpiifort

# etc...
```
