
![logo](images/logo_ndm.png)


## NDM

![image](images/cea_small.png)
[DES/ISAS/DMN/SRMP](https://www.universite-paris-saclay.fr/laboratoires/service-de-recherches-de-metallurgie-physique-des/isas/dmn)

### Overview

``NDM`` is fortran code computing empirical potential molecular dynamics (dynamique moléculaire en potentiels empiriques)
* etc TODO

See [NDM TODO better introduction](https://inis.iaea.org/collection/NCLCollectionStore/_Public/37/064/37064766.pdf)


NDM is built around
- [Fortran 2008+](https://fr.wikipedia.org/wiki/Fortran)
- [CMake](https://cmake.org/)
- [GNU compilers products](https://www.intel.com/content/www/us/en/developer/tools/oneapi/overview.html).
- [OneApi Intel products](https://www.intel.com/content/www/us/en/developer/tools/oneapi/overview.html).

All the input data are prescribed non-interactively,
using command lines and/or ASCII files.

NDM can be compiled and run on any Unix-like system (and TODO windows).


### Installation

```
git clone --branch ndm2021_cv ssh://gitolite@ssh-codev-tuleap.intra.cea.fr:2044/ndm/NDM.git NDM
```

### Compilation

### Compiling NDM with Cmake

Prerequisities
It is recommended to use Intel's oneAPI suite including Intel's mpi implementation, the mpiifort wrapper and MKL.
To use lammps it must be compiled as a shared library. The shared library Lammps compilation should create the required liblammps.so dynamic library. It is recommended to compile Lammps with the same compiler suite as NDM. see lammps doc

Standard NDM compilation using cmake
For convenience the following steps have been written in a bash script which you should edit according to your own settings.
1. Choose the relevant preset file with your prefered NDM configuration (using MPI/LAMMPS).
load the cache file:
cmake -C cache_file

2. then generate NDM cmake configuration (from any directory) :
cmake -B BUI_DIR -S SRC_DIR
with BUI_DIR the directory where you want NDM to be built and SRC_DIR the directory containing the main CMakeLists.txt file (usually NDM, parent from the src directory containing .F90 sources).

3. from the BUI_DIR, build NDM :
cmake --build .
the compilation can be parallelized using the --parallel <n_proc> option.

*. When Cmake fails to find the relevant dependencies it is often required to define additionnal variables to indiquate the libraries locations such as MPI_HOME, MPI_ROOT, etc... They are listed in the bash script.





Compilation instructions are provided in the documentation included in the
distribution repository
- file [README_ndm_compilation.md](READMES/README_ndm_compilation.md)

Also see TODO
- github documentation [ndm-docs TODO](https://jpc.github.io/ndm-docs/contents/installation.html).


### Integration Tests TODO

- Data are small,
  located at [NDM/examples directories](https://codev-tuleap.intra.cea.fr/plugins/git/ndm/NDM) branch `ndm2021_cv`.

- Data are big,
  located in a [NDM_TESTS separate git repository](https://codev-tuleap.intra.cea.fr/plugins/git/ndm/NDM_TESTS.git) TODO.
