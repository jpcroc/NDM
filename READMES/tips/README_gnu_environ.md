
### GNU compilation

On centos-7 (for example) user have to set *new* version 
of gnu compilers to compile MILADY.

#### Using `module` command

On clusters (usually) user could use `module`

```
module avail gcc
module load gcc/9.3.0

module avail openmpi
module load openmpi/gcc_9.3.0/4.0.1
etc...
```

#### Using `scl` command

On other local hosts, sometimes, user could use `scl`

```
# gnu 10 case
scl enable devtoolset-10 bash
```

```
# example on centos-7 is221713
which gfortran
  /usr/bin/gfortran

gfortran -v
  Utilisation des specs internes.
  COLLECT_GCC=gfortran
  COLLECT_LTO_WRAPPER=/usr/libexec/gcc/x86_64-redhat-linux/4.8.5/lto-wrapper
  Target: x86_64-redhat-linux
  Configuré avec: etc...
  Modèle de thread: posix
  gcc version 4.8.5 20150623 (Red Hat 4.8.5-44) (GCC) 


scl enable devtoolset-10 bash

which gfortran
  /opt/rh/devtoolset-10/root/usr/bin/gfortran

gfortran -v
  Using built-in specs.
  COLLECT_GCC=gfortran
  COLLECT_LTO_WRAPPER=/opt/rh/devtoolset-10/root/usr/libexec/gcc/x86_64-redhat-linux/10/lto-wrapper
  Target: x86_64-redhat-linux
  Configured with: etc...
  Thread model: posix
  Supported LTO compression algorithms: zlib
  gcc version 10.2.1 20200804 (Red Hat 10.2.1-2) (GCC) 
```
